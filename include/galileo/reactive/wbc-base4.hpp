#pragma once

#include "galileo/reactive/fwd.hpp"
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

        using scalar_t = double;
        using vector_t = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
        using matrix_t = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

        using Vector6 = Eigen::Matrix<scalar_t, 6, 1>;
        using Matrix6 = Eigen::Matrix<scalar_t, 6, 6>;

        using Model = pinocchio::ModelTpl<scalar_t, 0, pinocchio::JointCollectionDefaultTpl>;
        using Data = pinocchio::DataTpl<scalar_t, 0, pinocchio::JointCollectionDefaultTpl>;

        struct Info
        {
            int nq = 19;
            int nv = 18;
            int actuatedDofNum = 12;
            int numThreeDofContacts = 4;
        }; // struct Info

        class Task
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            Task() = default;

            Task(matrix_t a, vector_t b, matrix_t d, vector_t f) : a_(std::move(a)), d_(std::move(d)), b_(std::move(b)), f_(std::move(f)) {}

            explicit Task(size_t numDecisionVars)
                : Task(matrix_t::Zero(0, numDecisionVars), vector_t::Zero(0), matrix_t::Zero(0, numDecisionVars), vector_t::Zero(0)) {}

            Task operator+(const Task &rhs) const
            {
                return {concatenateMatrices(a_, rhs.a_), concatenateVectors(b_, rhs.b_), concatenateMatrices(d_, rhs.d_),
                        concatenateVectors(f_, rhs.f_)};
            }

            Task operator*(scalar_t rhs) const
            { // clang-format off
       return {a_.cols() > 0 ? rhs * a_ : a_,
               b_.cols() > 0 ? rhs * b_ : b_,
               d_.cols() > 0 ? rhs * d_ : d_,
               f_.cols() > 0 ? rhs * f_ : f_}; // clang-format on
            }

            matrix_t a_, d_;
            vector_t b_, f_;

            static matrix_t concatenateMatrices(matrix_t m1, matrix_t m2)
            {
                if (m1.cols() <= 0)
                {
                    return m2;
                }
                else if (m2.cols() <= 0)
                {
                    return m1;
                }
                assert(m1.cols() == m2.cols());
                matrix_t res(m1.rows() + m2.rows(), m1.cols());
                res << m1, m2;
                return res;
            }

            static vector_t concatenateVectors(const vector_t &v1, const vector_t &v2)
            {
                if (v1.cols() <= 0)
                {
                    return v2;
                }
                else if (v2.cols() <= 0)
                {
                    return v1;
                }
                assert(v1.cols() == v2.cols());
                vector_t res(v1.rows() + v2.rows());
                res << v1, v2;
                return res;
            }
        }; // class Task

        class WBCBase
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            WBCBase() {}

            WBCBase(Model model, const Info &info, const std::vector<multibody::EndEffector> &ees)
                : model_(model),
                  data_measured_(model_),
                  data_desired_(model),
                  info_(info),
                  ees_(ees)
            {
                numDecisionVars_ = info_.nv + 3 * info_.numThreeDofContacts + info_.actuatedDofNum;

                qMeasured_ = vector_t(info_.nq);
                vMeasured_ = vector_t(info_.nv);

                torqueLimits_ = vector_t::Constant(3, 24);

                std::cout << "WBCBase created" << std::endl;
            }

            vector_t update(const vector_t &stateDesired, const vector_t &inputDesired, const vector_t &rbdStateMeasured, size_t mode, scalar_t period)
            {
                auto start = std::chrono::high_resolution_clock::now();
                numContacts_ = info_.numThreeDofContacts; // for now
                contactFlag_ = {true, true, true, true};  // for now

                updateMeasured(rbdStateMeasured);
                updateDesired(stateDesired, inputDesired);

                // Constraints
                Task constraints = formulateConstraints();
                size_t numConstraints = constraints.b_.size() + constraints.f_.size();

                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A(numConstraints, getNumDecisionVars());
                vector_t lbA(numConstraints), ubA(numConstraints);
                A << constraints.a_,
                    constraints.d_;

                lbA << constraints.b_,
                    -qpOASES::INFTY * vector_t::Ones(constraints.f_.size());
                ubA << constraints.b_,
                    constraints.f_;

                // Cost
                Task weighedTask = formulateWeightedTasks(stateDesired, inputDesired, period);
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> H = weighedTask.a_.transpose() * weighedTask.a_;
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
                auto qpProblem = qpOASES::QProblem(getNumDecisionVars(), numConstraints);
                qpOASES::Options options;
                options.setToMPC();
                options.printLevel = qpOASES::PL_LOW;
                options.enableEqualities = qpOASES::BT_TRUE;
                qpProblem.setOptions(options);
                int nWsr = 20;

                qpProblem.init(H.data(), g.data(), A.data(), nullptr, nullptr, lbA.data(), ubA.data(), nWsr);
                vector_t qpSol(getNumDecisionVars());

                qpProblem.getPrimalSolution(qpSol.data());

                end = std::chrono::high_resolution_clock::now();
                std::cout << "QP Solve Time Taken: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microseconds" << std::endl;

                return qpSol;
            }

            void updateMeasured(const vector_t &rbdStateMeasured)
            {
                qMeasured_ = rbdStateMeasured.head(info_.nq);
                vMeasured_ = rbdStateMeasured.tail(info_.nv);

                // For floating base EoM task
                pinocchio::forwardKinematics(model_, data_measured_, qMeasured_, vMeasured_);
                pinocchio::computeJointJacobians(model_, data_measured_);
                pinocchio::updateFramePlacements(model_, data_measured_);
                pinocchio::crba(model_, data_measured_, qMeasured_);
                data_measured_.M.triangularView<Eigen::StrictlyLower>() = data_measured_.M.transpose().triangularView<Eigen::StrictlyLower>();
                pinocchio::nonLinearEffects(model_, data_measured_, qMeasured_, vMeasured_);
                j_ = matrix_t(3 * info_.numThreeDofContacts, info_.nv);
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    Eigen::Matrix<scalar_t, 6, Eigen::Dynamic> jac;
                    jac.setZero(6, info_.nv);
                    pinocchio::getFrameJacobian(model_, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
                    // std::cout << "jac for " << ees_[i].frame_name << ":\n" << jac << std::endl;
                    j_.block(3 * i, 0, 3, info_.nv) = jac.template topRows<3>();
                }

                // For not contact motion task
                pinocchio::computeJointJacobiansTimeVariation(model_, data_measured_, qMeasured_, vMeasured_);
                dj_ = matrix_t(3 * info_.numThreeDofContacts, info_.nv);
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    Eigen::Matrix<scalar_t, 6, Eigen::Dynamic> jac;
                    jac.setZero(6, info_.nv);
                    pinocchio::getFrameJacobianTimeVariation(model_, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED, jac);
                    // std::cout << "djac for " << ees_[i].frame_name << ":\n" << jac << std::endl;
                    dj_.block(3 * i, 0, 3, info_.nv) = jac.template topRows<3>();
                }
            }

            void updateDesired(const vector_t &stateDesired, const vector_t &inputDesired)
            {
                const vector_t qDesired = stateDesired.head(info_.nq); // some mapping
                pinocchio::forwardKinematics(model_, data_desired_, qDesired);
                pinocchio::computeJointJacobians(model_, data_desired_, qDesired);
                pinocchio::updateFramePlacements(model_, data_desired_);
                const vector_t vDesired = stateDesired.tail(info_.nv); // some mapping
                pinocchio::forwardKinematics(model_, data_desired_, qDesired, vDesired);

                pinocchio::centerOfMass(model_, data_desired_, qDesired, false);
                pinocchio::computeCentroidalMap(model_, data_desired_, qDesired);
            }

            Task formulateConstraints()
            {
                return formulateFloatingBaseEomTask() + formulateTorqueLimitsTask() + formulateFrictionConeTask() + formulateNoContactMotionTask();
                // return formulateFloatingBaseEomTask() + formulateFrictionConeTask() + formulateNoContactMotionTask();
            }

            Task formulateWeightedTasks(const vector_t &stateDesired, const vector_t &inputDesired, scalar_t period)
            {
                return formulateSwingLegTask() * weightSwingLeg_ + formulateBaseAccelTask(stateDesired, inputDesired, period) * weightBaseAccel_ +
                       formulateContactForceTask(inputDesired) * weightContactForce_;
                // return formulateBaseAccelTask(stateDesired, inputDesired, period) * weightBaseAccel_ + formulateContactForceTask(inputDesired) * weightContactForce_;
                // return formulateTrackingTask(stateDesired, inputDesired, period);
            }

            // a: [nv, numDecisionVars_]
            // b: [nv]
            Task formulateFloatingBaseEomTask()
            {

                matrix_t s(info_.actuatedDofNum, info_.nv);
                s.block(0, 0, info_.actuatedDofNum, 6).setZero();
                s.block(0, 6, info_.actuatedDofNum, info_.actuatedDofNum).setIdentity();

                matrix_t a = (matrix_t(info_.nv, numDecisionVars_) << data_measured_.M, -j_.transpose(), -s.transpose()).finished();
                vector_t b = -data_measured_.nle;

                return {a, b, matrix_t(), vector_t()};
            }

            // d: [2 * actuatedDofNum, numDecisionVars_]
            // f: [2 * actuatedDofNum]
            Task formulateTorqueLimitsTask()
            {
                matrix_t d(2 * info_.actuatedDofNum, numDecisionVars_);
                d.setZero();
                matrix_t i = matrix_t::Identity(info_.actuatedDofNum, info_.actuatedDofNum);
                d.block(0, info_.nv + 3 * info_.numThreeDofContacts, info_.actuatedDofNum, info_.actuatedDofNum) = i;
                d.block(info_.actuatedDofNum, info_.nv + 3 * info_.numThreeDofContacts, info_.actuatedDofNum,
                        info_.actuatedDofNum) = -i;
                vector_t f(2 * info_.actuatedDofNum);
                for (size_t l = 0; l < 2 * info_.actuatedDofNum / 3; ++l)
                {
                    f.segment<3>(3 * l) = torqueLimits_;
                }

                return {matrix_t(), vector_t(), d, f};
            }

            // a: [3 * numContacts_, numDecisionVars_]
            // b: [3 * numContacts_]
            Task formulateNoContactMotionTask()
            {
                matrix_t a(3 * numContacts_, numDecisionVars_);
                vector_t b(a.rows());
                a.setZero();
                b.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; i++)
                {
                    // if (contactFlag_[i])
                    // {
                        a.block(3 * j, 0, 3, info_.nv) = j_.block(3 * i, 0, 3, info_.nv);
                        b.segment(3 * j, 3) = -dj_.block(3 * i, 0, 3, info_.nv) * vMeasured_;
                        j++;
                    // }
                }

                return {a, b, matrix_t(), vector_t()};
            }

            // a: [3 * (numThreeDofContacts - numContacts_), numDecisionVars_]
            // b: [3 * (numThreeDofContacts - numContacts_)]
            // d: [5 * numContacts_ + 3 * (numThreeDofContacts - numContacts_), numDecisionVars_]
            // f: [5 * numContacts_ + 3 * (numThreeDofContacts - numContacts_)]
            Task formulateFrictionConeTask()
            {
                matrix_t a(3 * (info_.numThreeDofContacts - numContacts_), numDecisionVars_);
                a.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    if (!contactFlag_[i])
                    {
                        a.block(3 * j++, info_.nv + 3 * i, 3, 3) = matrix_t::Identity(3, 3);
                    }
                }
                vector_t b(a.rows());
                b.setZero();

                matrix_t frictionPyramic(5, 3);
                frictionPyramic << 0, 0, -1,
                    1, 0, -frictionCoeff_,
                    -1, 0, -frictionCoeff_,
                    0, 1, -frictionCoeff_,
                    0, -1, -frictionCoeff_;

                matrix_t d(5 * numContacts_ + 3 * (info_.numThreeDofContacts - numContacts_), numDecisionVars_);
                d.setZero();
                j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    if (contactFlag_[i])
                    {
                        d.block(5 * j++, info_.nv + 3 * i, 5, 3) = frictionPyramic;
                    }
                }
                vector_t f = Eigen::VectorXd::Zero(d.rows());

                return {a, b, d, f};
            }

            // a: [numDecisionVars_, numDecisionVars_]
            // b: [numDecisionVars_]
            Task formulateTrackingTask(const vector_t &stateDesired, const vector_t &inputDesired, scalar_t period)
            {
                // identity tracking task
                matrix_t a = matrix_t::Identity(numDecisionVars_, numDecisionVars_);
                vector_t b = inputDesired;

                return {a, b, matrix_t(), vector_t()};
            }

            // a: [6, numDecisionVars_]
            // b: [6]
            Task formulateBaseAccelTask(const vector_t &stateDesired, const vector_t &inputDesired, scalar_t period)
            {
                matrix_t a(6, numDecisionVars_);
                a.setZero();
                a.block(0, 0, 6, 6) = matrix_t::Identity(6, 6);

                vector_t jointAccel = inputDesired.head(info_.nv);

                const vector_t qDesired = stateDesired.head(info_.nq);
                const vector_t vDesired = stateDesired.tail(info_.nv);

                const auto &A = data_desired_.Ag;
                const Matrix6 Ab = A.template leftCols<6>();
                const auto AbInv = computeFloatingBaseCentroidalMomentumMatrixInverse(Ab);
                const auto Aj = A.rightCols(info_.actuatedDofNum);
                const auto ADot = pinocchio::dccrba(model_, data_desired_, qDesired, vDesired);
                Vector6 centroidalMomentumRate = data_desired_.mass[0] * getNormalizedCentroidalMomentumRate(inputDesired);
                centroidalMomentumRate.noalias() -= ADot * vDesired;
                centroidalMomentumRate.noalias() -= Aj * jointAccel;

                Vector6 b = AbInv * centroidalMomentumRate;

                return {a, b, matrix_t(), vector_t()};
            }

            // a: [3 * (numThreeDofContacts - numContacts_), numDecisionVars_]
            // b: [3 * (numThreeDofContacts - numContacts_)]
            Task formulateSwingLegTask()
            {
                matrix_t a(3 * (info_.numThreeDofContacts - numContacts_), numDecisionVars_);
                vector_t b(a.rows());
                a.setZero();
                b.setZero();
                size_t j = 0;
                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    if (!contactFlag_[i])
                    {
                        vector_t pos_measured_i = data_measured_.oMf[ees_[i].frame_idx].translation();
                        vector_t vel_measured_i = pinocchio::getFrameVelocity(model_, data_measured_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED).linear();

                        vector_t pos_desired_i = data_desired_.oMf[ees_[i].frame_idx].translation();
                        vector_t vel_desired_i = pinocchio::getFrameVelocity(model_, data_desired_, ees_[i].frame_idx, pinocchio::LOCAL_WORLD_ALIGNED).linear();

                        vector_t accel = swingKp_ * (pos_desired_i - pos_measured_i) + swingKd_ * (vel_desired_i - vel_measured_i);
                        a.block(3 * j, 0, 3, info_.nv) = j_.block(3 * i, 0, 3, info_.nv);
                        b.segment(3 * j, 3) = accel - dj_.block(3 * i, 0, 3, info_.nv) * vMeasured_;
                        j++;
                    }
                }

                return {a, b, matrix_t(), vector_t()};
            }

            // a: [3 * numContacts_, numDecisionVars_]
            // b: [3 * numContacts_]
            Task formulateContactForceTask(const vector_t &inputDesired) const
            {
                matrix_t a(3 * info_.numThreeDofContacts, numDecisionVars_);
                vector_t b(a.rows());
                a.setZero();

                for (size_t i = 0; i < info_.numThreeDofContacts; ++i)
                {
                    a.block(3 * i, info_.nv + 3 * i, 3, 3) = matrix_t::Identity(3, 3);
                }
                b = inputDesired.head(a.rows());

                return {a, b, matrix_t(), vector_t()};
            }

            template <typename SCALAR_T>
            Eigen::Matrix<SCALAR_T, 6, 1> getNormalizedCentroidalMomentumRate(const Eigen::Matrix<SCALAR_T, Eigen::Dynamic, 1> &input)
            {
                const Eigen::Matrix<SCALAR_T, 3, 1> gravityVector(SCALAR_T(0.0), SCALAR_T(0.0), SCALAR_T(-9.81));
                Eigen::Matrix<SCALAR_T, 6, 1> centroidalMomentumRate;
                centroidalMomentumRate << data_desired_.mass[0] * gravityVector, Eigen::Matrix<SCALAR_T, 3, 1>::Zero();

                for (size_t i = 0; i < info_.numThreeDofContacts; i++)
                {
                    const auto contactForceInWorldFrame = input.template segment<3>(info_.nv + 3 * i);
                    const auto positionComToContactPointInWorldFrame = (data_desired_.oMf[ees_[i].frame_idx].translation() - data_desired_.com[0]);
                    centroidalMomentumRate.template head<3>() += contactForceInWorldFrame;
                    centroidalMomentumRate.template tail<3>().noalias() += positionComToContactPointInWorldFrame.cross(contactForceInWorldFrame);
                } // end of i loop

                // normalize by the total mass
                centroidalMomentumRate /= data_desired_.mass[0];

                return centroidalMomentumRate;
            }

            template <typename SCALAR_T>
            Eigen::Matrix<SCALAR_T, 6, 6> computeFloatingBaseCentroidalMomentumMatrixInverse(const Eigen::Matrix<SCALAR_T, 6, 6> &Ab)
            {
                const SCALAR_T mass = Ab(0, 0);
                Eigen::Matrix<SCALAR_T, 3, 3> Ab_22_inv = Ab.template block<3, 3>(3, 3).inverse();
                Eigen::Matrix<SCALAR_T, 6, 6> Ab_inv = Eigen::Matrix<SCALAR_T, 6, 6>::Zero();
                Ab_inv << 1.0 / mass * Eigen::Matrix<SCALAR_T, 3, 3>::Identity(), -1.0 / mass * Ab.template block<3, 3>(0, 3) * Ab_22_inv,
                    Eigen::Matrix<SCALAR_T, 3, 3>::Zero(), Ab_22_inv;
                return Ab_inv;
            }

            size_t getNumDecisionVars() const { return numDecisionVars_; }

        protected:
            Model model_;

            Data data_measured_;
            Data data_desired_;

            matrix_t j_;
            matrix_t dj_;

            vector_t qMeasured_;
            vector_t vMeasured_;

            vector_t torqueLimits_;

            scalar_t frictionCoeff_ = 0.4;
            scalar_t swingKp_ = 350.;
            scalar_t swingKd_ = 37.;

            scalar_t weightSwingLeg_ = 100.;
            scalar_t weightBaseAccel_ = 1.;
            scalar_t weightContactForce_ = 0.01;

            Info info_;

            std::vector<multibody::EndEffector> ees_;

            std::array<bool, 4> contactFlag_{};

            size_t numDecisionVars_;
            size_t numContacts_{};

        }; // class WBCBase

    } // namespace reactive
} // namespace galileo