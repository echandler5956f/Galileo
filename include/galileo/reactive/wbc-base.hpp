#ifndef __galileo_reactive_wbc_base_hpp__
#define __galileo_reactive_wbc_base_hpp__

#include "galileo/reactive/fwd.hpp"
#include "galileo/reactive/task.hpp"
#include "galileo/multibody/end-effectors.hpp"

#include <pinocchio/fwd.hpp>
#include <pinocchio/algorithm/fwd.hpp>
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/parsers/mjcf.hpp>

#include <qpOASES.hpp>
#include <memory>
#include <iostream>
#include <chrono>

namespace galileo
{
    namespace reactive
    {
        struct Info
        {
            int nq;
            int nv;
            int actuatedDofNum;
            int numThreeDofContacts;
        }; // struct Info

        // decision vars are [a, f, tau]
        template <typename NumScalar, int Options = 0>
        class WBCBase
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using Model = pinocchio::ModelTpl<NumScalar, 0, pinocchio::JointCollectionDefaultTpl>;
            using Data = pinocchio::DataTpl<NumScalar, 0, pinocchio::JointCollectionDefaultTpl>;

            using SE3_t = pinocchio::SE3Tpl<NumScalar, Options>;

            using vector3_t = Eigen::Matrix<NumScalar, 3, 1>;
            using vector6_t = Eigen::Matrix<NumScalar, 6, 1>;
            using matrix6_t = Eigen::Matrix<NumScalar, 6, 6>;
            using matrix_t = Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic, Options>;
            using vector_t = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1, Options>;

            WBCBase()
            {
            }

            WBCBase(Model model, const Info &info, const std::vector<multibody::EndEffector> &ees)
                : model_(model),
                  data_measured_(model_),
                  data_desired_(model_),
                  info_(info),
                  ees_(ees)
            {
                num_decision_vars_ = info_.nv + 3 * info_.numThreeDofContacts + info_.actuatedDofNum;

                num_contacts_ = info_.numThreeDofContacts; // for now
                // all contacts are active
                contact_flag_ = std::vector<bool>(info_.numThreeDofContacts, true); // for now

                q_measured_ = vector_t(info_.nq);
                v_measured_ = vector_t(info_.nv);

                control_limits_ = vector3_t::Constant(3, 24);

                std::cout << "WBCBase created" << std::endl;
            }

            template <typename DesStateVector, typename DesControlVector, typename RBDStateVector>
            vector_t update(const Eigen::MatrixBase<DesStateVector> &stateDesired, const Eigen::MatrixBase<DesControlVector> &controlDesired, const Eigen::MatrixBase<RBDStateVector> &rbdStateMeasured, size_t mode, NumScalar period)
            {
                auto start = std::chrono::high_resolution_clock::now();
                updateMeasured(rbdStateMeasured.derived());
                updateDesired(stateDesired.derived(), controlDesired.derived());

                TaskDefault<NumScalar> constraints = std::move(formulateConstraints());
                // TaskDefault<NumScalar> constraints = std::move(formulateFloatingBaseEomTask());

                size_t numConstraints = constraints.b_.size() + constraints.f_.size();

                Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(numConstraints, getNumDecisionVars());
                vector_t lbA(numConstraints), ubA(numConstraints);
                A << constraints.a_,
                    constraints.d_;

                lbA << constraints.b_,
                    -qpOASES::INFTY * vector_t::Ones(constraints.f_.size());
                ubA << constraints.b_,
                    constraints.f_;

                // Cost
                TaskDefault<NumScalar> weighedTask = std::move(formulateWeightedTasks(stateDesired, controlDesired, period));
                Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> H = weighedTask.a_.transpose() * weighedTask.a_;
                vector_t g = -weighedTask.a_.transpose() * weighedTask.b_;

                // std::cout << "A: " << A << std::endl;
                // std::cout << "lbA: " << lbA << std::endl;
                // std::cout << "ubA: " << ubA << std::endl;
                // std::cout << "H: " << H << std::endl;
                // std::cout << "g: " << g << std::endl;

                auto end = std::chrono::high_resolution_clock::now();
                std::cout << "QP Setup Time Taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;

                start = std::chrono::high_resolution_clock::now();

                // Solve
                qpOASES::QProblem qpProblem = qpOASES::QProblem(getNumDecisionVars(), numConstraints);
                qpOASES::Options options;
                options.setToMPC();
                options.printLevel = qpOASES::PL_LOW;
                options.enableEqualities = qpOASES::BT_TRUE;
                qpProblem.setOptions(options);
                int nWsr = 20;

                // std::cout << "Solving QP" << std::endl;

                qpProblem.init(H.data(), g.data(), A.data(), nullptr, nullptr, lbA.data(), ubA.data(), nWsr);
                vector_t qpSol(getNumDecisionVars());

                qpProblem.getPrimalSolution(qpSol.data());

                // std::cout << "QP solved" << std::endl;

                // vector_t tau = qpSol.tail(info_.actuatedDofNum);

                end = std::chrono::high_resolution_clock::now();
                std::cout << "QP Solve Time Taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;

                // return tau;
                return qpSol;
            }

        protected:
            template <typename RBDStateVector>
            void updateMeasured(const Eigen::MatrixBase<RBDStateVector> &rbdStateMeasured)
            {
                q_measured_ = rbdStateMeasured.head(info_.nq);
                v_measured_ = rbdStateMeasured.tail(info_.nv);

                Model model = model_;

                // For floating base EoM task
                pinocchio::forwardKinematics(model, data_measured_, q_measured_, v_measured_);
                pinocchio::computeJointJacobians(model, data_measured_, q_measured_);
                pinocchio::updateFramePlacements(model, data_measured_);
                pinocchio::crba(model, data_measured_, q_measured_);
                data_measured_.M.template triangularView<Eigen::StrictlyLower>() = data_measured_.M.transpose().template triangularView<Eigen::StrictlyLower>();
                pinocchio::nonLinearEffects(model, data_measured_, q_measured_, v_measured_);
                frame_jac_ = matrix_t(3 * info_.numThreeDofContacts, info_.nv);
                frame_jac_.setZero();
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    matrix_t jac = matrix_t::Zero(6, info_.nv);
                    pinocchio::getFrameJacobian(model, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
                    frame_jac_.block(3 * i, 0, 3, info_.nv) = jac.template topRows<3>();
                }

                // For not contact motion task
                pinocchio::computeJointJacobiansTimeVariation(model, data_measured_, q_measured_, v_measured_);
                frame_jac_dot_ = matrix_t(3 * info_.numThreeDofContacts, info_.nv);
                frame_jac_dot_.setZero();
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    matrix_t jac = matrix_t::Zero(6, info_.nv);
                    pinocchio::getFrameJacobianTimeVariation(model, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
                    frame_jac_dot_.block(3 * i, 0, 3, info_.nv) = jac.template topRows<3>();
                }

                body_jac_ = matrix_t(6, info_.nv);
                body_jac_.setZero();
                pinocchio::getFrameJacobian(model, data_measured_, model.getFrameId("base", pinocchio::BODY), pinocchio::LOCAL_WORLD_ALIGNED, body_jac_);
                body_jac_dot_ = matrix_t(6, info_.nv);
                body_jac_dot_.setZero();
                pinocchio::getFrameJacobianTimeVariation(model, data_measured_, model.getFrameId("base", pinocchio::BODY), pinocchio::LOCAL_WORLD_ALIGNED, body_jac_dot_);
            }

            template <typename DesStateVector, typename DesControlVector>
            void updateDesired(const Eigen::MatrixBase<DesStateVector> &stateDesired, const Eigen::MatrixBase<DesControlVector> &controlDesired)
            {
                Model model = model_;

                const vector_t qDesired = stateDesired.head(info_.nq); // some mapping

                pinocchio::forwardKinematics(model, data_desired_, qDesired);
                pinocchio::computeJointJacobians(model, data_desired_, qDesired);
                pinocchio::updateFramePlacements(model, data_desired_);
                // updateCentroidalDynamics(pinocchioInterfaceDesired_, info_, qDesired);
                const vector_t vDesired = stateDesired.tail(info_.nv); // some mapping
                pinocchio::forwardKinematics(model, data_desired_, qDesired, vDesired);
            }

            TaskDefault<NumScalar> formulateConstraints()
            {
                return taskVertcat(
                    std::move(formulateFloatingBaseEomTask()),
                    std::move(taskVertcat(
                        std::move(formulateTorqueLimitsTask()),
                        std::move(taskVertcat(
                            std::move(formulateFrictionConeTask()),
                            std::move(formulateNoContactMotionTask()))))));
                // return taskVertcat(
                //     std::move(formulateFloatingBaseEomTask()),
                //     std::move(taskVertcat(
                //         std::move(formulateFrictionConeTask()),
                //         std::move(formulateNoContactMotionTask()))));
                // return taskVertcat(
                //     std::move(formulateFloatingBaseEomTask()),
                //     std::move(formulateFrictionConeTask()));
                    // return taskVertcat(
                    //     std::move(formulateFloatingBaseEomTask()),
                    //     std::move(formulateNoContactMotionTask()));

                // return std::move(formulateFloatingBaseEomTask());
            }

            template <typename DesStateVector, typename DesControlVector>
            TaskDefault<NumScalar> formulateWeightedTasks(const Eigen::MatrixBase<DesStateVector> &stateDesired, const Eigen::MatrixBase<DesControlVector> &controlDesired, NumScalar period)
            {
                return taskVertcat(
                    std::move(formulateSwingLegTask() * weightSwingLeg_),
                    std::move(taskVertcat(
                        std::move(formulateBaseAccelTask(stateDesired.derived(), controlDesired.derived(), period) * weightBaseAccel_),
                        std::move(formulateContactForceTask(controlDesired.derived()) * weightContactForce_))));
                // return taskVertcat(
                //     std::move(formulateBaseAccelTask(stateDesired.derived(), controlDesired.derived(), period) * weightBaseAccel_),
                //     std::move(formulateContactForceTask(controlDesired.derived()) * weightContactForce_));
                // return std::move(formulateContactForceTask(controlDesired.derived()) * weightContactForce_);
                // return std::move(formulateBaseAccelTask(stateDesired.derived(), controlDesired.derived(), period) * weightBaseAccel_);
            }

            TaskDefault<NumScalar> formulateFloatingBaseEomTask()
            {
                matrix_t s(info_.actuatedDofNum, info_.nv);
                s.setZero();
                s.block(0, 6, info_.actuatedDofNum, info_.actuatedDofNum).setIdentity();

                // matrix_t a(info_.nv, num_decision_vars_);

                // a.setZero();
                // a.block(0, 0, info_.nv, info_.nv) = data_measured_.M;
                // a.block(0, info_.nv, info_.nv, 3 * info_.numThreeDofContacts) = -frame_jac_.transpose();
                // a.block(0, info_.nv + 3 * info_.numThreeDofContacts, info_.nv, info_.actuatedDofNum) = -s.transpose();

                // should be equal to this
                matrix_t a = (matrix_t(info_.nv, num_decision_vars_) << data_measured_.M, -frame_jac_.transpose(), -s.transpose()).finished();

                vector_t b = -data_measured_.nle;

                return TaskDefault<NumScalar>(a, b, matrix_t(), vector_t());
            }

            TaskDefault<NumScalar> formulateTorqueLimitsTask()
            {
                matrix_t d(2 * info_.actuatedDofNum, num_decision_vars_);
                d.setZero();
                matrix_t i = matrix_t::Identity(info_.actuatedDofNum, info_.actuatedDofNum);
                d.block(0, info_.nv + 3 * info_.numThreeDofContacts, info_.actuatedDofNum, info_.actuatedDofNum) = i;
                d.block(info_.actuatedDofNum, info_.nv + 3 * info_.numThreeDofContacts, info_.actuatedDofNum,
                        info_.actuatedDofNum) = -i;
                vector_t f(2 * info_.actuatedDofNum);
                f.setZero();
                for (size_t l = 0; l < 2 * info_.actuatedDofNum / 3; ++l)
                {
                    f.segment(3 * l, 3) = control_limits_;
                }

                return TaskDefault<NumScalar>(matrix_t(), vector_t(), d, f);
            }

            TaskDefault<NumScalar> formulateNoContactMotionTask()
            {
                matrix_t a(3 * num_contacts_, num_decision_vars_);
                vector_t b(a.rows());
                a.setZero();
                b.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; i++)
                {
                    if (contact_flag_[i])
                    {
                        a.block(3 * j, 0, 3, info_.nv) = frame_jac_.block(3 * i, 0, 3, info_.nv);
                        b.segment(3 * j, 3) = -frame_jac_dot_.block(3 * i, 0, 3, info_.nv) * v_measured_;
                        j++;
                    }
                }

                return TaskDefault<NumScalar>(a, b, matrix_t(), vector_t());
            }

            TaskDefault<NumScalar> formulateFrictionConeTask()
            {
                matrix_t a(3 * (info_.numThreeDofContacts - num_contacts_), num_decision_vars_);
                a.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    if (!contact_flag_[i])
                    {
                        a.block(3 * j++, info_.nv + 3 * i, 3, 3) = matrix_t::Identity(3, 3);
                    }
                }
                vector_t b(a.rows());
                b.setZero();

                matrix_t frictionPyramic(5, 3);
                frictionPyramic << 0, 0, -1,
                    1, 0, -friction_coeff_,
                    -1, 0, -friction_coeff_,
                    0, 1, -friction_coeff_,
                    0, -1, -friction_coeff_;

                matrix_t d(5 * num_contacts_ + 3 * (info_.numThreeDofContacts - num_contacts_), num_decision_vars_);
                d.setZero();
                j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    if (contact_flag_[i])
                    {
                        d.block(5 * j++, info_.nv + 3 * i, 5, 3) = frictionPyramic;
                    }
                }
                vector_t f = Eigen::VectorXd::Zero(d.rows());

                return TaskDefault<NumScalar>(a, b, d, f);
            }

            TaskDefault<NumScalar> formulateBaseAccelTask(const vector_t &stateDesired, const vector_t &controlDesired, NumScalar period)
            {
                Model model = model_;

                matrix_t a(6, num_decision_vars_);
                a.setZero();
                a.block(0, 0, 6, info_.nv) = body_jac_;

                vector6_t b = vector6_t::Zero(6);
                vector3_t pos_desired = stateDesired.head(3);
                vector6_t vel_desired = stateDesired.segment(info_.nq, 6);

                // quaternion expects w, x, y, z, but stores internally as x, y, z, w (why?!)
                // SE3_t desired_frame = SE3_t(typename SE3_t::Quaternion(stateDesired(6), stateDesired(3), stateDesired(4), stateDesired(5)), pos_desired);
                // SE3_t iMd = data_measured_.oMf[0].actInv(desired_frame);
                // vector6_t pose_err = pinocchio::log6(iMd).toVector();

                vector3_t pos_error = pos_desired - data_measured_.oMf[0].translation();
                vector6_t pose_err(pos_error(0), pos_error(1), pos_error(2), 0, 0, 0); // simplified for now

                vector_t vel_measured = pinocchio::getFrameVelocity(model, data_measured_, model.getFrameId("base", pinocchio::BODY), pinocchio::LOCAL_WORLD_ALIGNED).toVector();
                vector_t accel = body_kp_ * pose_err + body_kd_ * (vel_desired - vel_measured);
                b = accel - body_jac_dot_ * v_measured_;

                // std::cout << "a: \n" << a << std::endl;
                // std::cout << "b: \n" << b << std::endl;

                return TaskDefault<NumScalar>(a, b, matrix_t(), vector_t());
            }

            TaskDefault<NumScalar> formulateSwingLegTask()
            {
                Model model = model_;

                matrix_t a(3 * (info_.numThreeDofContacts - num_contacts_), num_decision_vars_);
                vector_t b(a.rows());
                a.setZero();
                b.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    vector3_t pos_measured_i = data_measured_.oMf[ees_[i].frame_idx].translation();
                    vector3_t vel_measured_i = pinocchio::getFrameVelocity(model, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED).linear();

                    vector3_t pos_desired_i = data_desired_.oMf[ees_[i].frame_idx].translation();
                    vector3_t vel_desired_i = pinocchio::getFrameVelocity(model, data_desired_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED).linear();
                    if (!contact_flag_[i])
                    {
                        vector3_t accel = swing_kp_ * (pos_desired_i - pos_measured_i) + swing_kd_ * (vel_desired_i - vel_measured_i);
                        a.block(3 * j, 0, 3, info_.nv) = frame_jac_.block(3 * i, 0, 3, info_.nv);
                        b.segment(3 * j, 3) = accel - frame_jac_dot_.block(3 * i, 0, 3, info_.nv) * v_measured_;
                        j++;
                    }
                }

                return TaskDefault<NumScalar>(a, b, matrix_t(), vector_t());
            }

            TaskDefault<NumScalar> formulateContactForceTask(const vector_t &controlDesired) const
            {
                // std::cout << "Formulating contact force task" << std::endl;
                matrix_t a(3 * info_.numThreeDofContacts, num_decision_vars_);
                vector_t b(a.rows());
                a.setZero();

                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    a.block(3 * i, info_.nv + 3 * i, 3, 3) = matrix_t::Identity(3, 3);
                }
                b = controlDesired.segment(info_.nv, 3 * info_.numThreeDofContacts);

                return TaskDefault<NumScalar>(a, b, matrix_t(), vector_t());
            }

            size_t getNumDecisionVars() const { return num_decision_vars_; }

            Model model_;

            Data data_measured_;
            Data data_desired_;

            matrix_t body_jac_;
            matrix_t body_jac_dot_;

            matrix_t frame_jac_;
            matrix_t frame_jac_dot_;

            vector_t q_measured_;
            vector_t v_measured_;

            vector3_t control_limits_;

            NumScalar friction_coeff_ = 0.6;
            NumScalar body_kp_ = 100.;
            NumScalar body_kd_ = 20.;
            NumScalar swing_kp_ = 350.;
            NumScalar swing_kd_ = 37.;

            NumScalar weightSwingLeg_ = 100.;
            NumScalar weightBaseAccel_ = 1.;
            NumScalar weightContactForce_ = 0.01;

            Info info_;

            std::vector<multibody::EndEffector> ees_;

            std::vector<bool> contact_flag_;

            size_t num_decision_vars_;
            size_t num_contacts_;

        }; // class WBCBase

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_wbc_base_hpp__