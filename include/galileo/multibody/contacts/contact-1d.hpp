#ifndef __galileo_multibody_contacts_contact_1d_hpp__
#define __galileo_multibody_contacts_contact_1d_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename PhaseSpec>
        struct Contact1dTpl;

        template <typename PhaseSpec>
        struct traits<Contact1dTpl<PhaseSpec>>
        {
            using PS = PhaseSpec;

            using ContactDataDerived = ContactData1dTpl<PhaseSpec>;
            using ContactModelDerived = ContactModel1dTpl<PhaseSpec>;

            static constexpr int NC = 1;

            // using RobotData_t = // TODO: add robot data
        };

        template <typename PhaseSpec>
        struct traits<ContactData1dTpl<PhaseSpec>>
        {
            using ContactDerived = Contact1dTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct traits<ContactModel1dTpl<PhaseSpec>>
        {
            using ContactDerived = Contact1dTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct ContactData1dTpl : public ContactDataBase<ContactData1dTpl<PhaseSpec>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using RobotData_t = typename PS::RobotData_t;
            using MatrixNcNv_t = typename PS::MatrixNcNv_t;
            using MatrixNcNdx_t = typename PS::MatrixNcNdx_t;
            using MatrixNcNu_t = typename PS::MatrixNcNu_t;

            RobotData_t *robot_data;
            pinocchio::FrameIndex frame;
            pinocchio::ReferenceFrame type;
            SE3 jMf;
            MatrixNcNv_t Jc;
            Force f;
            Force fext;
            MatrixNcNdx_t df_dx;
            MatrixNcNu_t df_du;

            typename SE3::ActionMatrixType fXj;
            VectorNc_t a0;
            MatrixNcNdx_t da0_dx;
            MatrixNv_t dtau_dq;

            // Need to define accessors for the base class

            // FORWARD_ACCESSOR(robot_data);
            // FORWARD_ACCESSOR(frame);
            // FORWARD_ACCESSOR(type);
            // FORWARD_ACCESSOR(jMf);
            // FORWARD_ACCESSOR(Jc);
            // FORWARD_ACCESSOR(f);
            // FORWARD_ACCESSOR(fext);
            // FORWARD_ACCESSOR(df_dx);
            // FORWARD_ACCESSOR(df_du);

            // FORWARD_ACCESSOR(fXj);
            // FORWARD_ACCESSOR(a0);
            // FORWARD_ACCESSOR(da0_dx);
            // FORWARD_ACCESSOR(dtau_dq);
        };

        template <typename PhaseSpec>
        struct ContactModel1dTpl : public ContactModelBase<ContactModel1dTpl<PhaseSpec>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactDerived = Contact1dTpl<PhaseSpec>;
            using ContactModelDerived = traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = traits<ContactDerived>::ContactDataDerived;

            static constexpr int NC = traits<ContactDerived>::NC;

            using MatrixNcNv_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NV, PS::Options>;
            using MatrixNcNdx_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NDX, PS::Options>;
            using MatrixNcNu_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NU, PS::Options>;

            template <typename StateVectorType>
            void calc(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                pinocchio::updateFramePlacement(*state_->get_pinocchio().get(), *data.pinocchio,
                                                id_);
                pinocchio::getFrameJacobian(*state_->get_pinocchio().get(), *data.pinocchio,
                                            id_, pinocchio::LOCAL, data.fJf);
                data.v = pinocchio::getFrameVelocity(*state_->get_pinocchio().get(),
                                                     *data.pinocchio, id_);

                data.a0_local =
                    pinocchio::getFrameClassicalAcceleration(
                        *state_->get_pinocchio().get(), *data.pinocchio, id_, pinocchio::LOCAL)
                        .linear();

                const Eigen::Ref<const Matrix3s> oRf = data.pinocchio->oMf[id_].rotation();
                if (gains_[0] != 0.)
                {
                    data.dp = data.pinocchio->oMf[id_].translation() -
                              (xref_ * Raxis_ * Vector3s::UnitZ());
                    data.dp_local.noalias() = oRf.transpose() * data.dp;
                    data.a0_local += gains_[0] * data.dp_local;
                }
                if (gains_[1] != 0.)
                {
                    data.a0_local += gains_[1] * data.v.linear();
                }
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    data.Jc.row(0) = (Raxis_ * data.fJf.template topRows<3>()).row(2);
                    data.a0[0] = (Raxis_ * data.a0_local)[2];
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    data.Jc.row(0) = (Raxis_ * oRf * data.fJf.template topRows<3>()).row(2);
                    data.a0[0] = (Raxis_ * oRf * data.a0_local)[2];
                    break;
                }
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                const pinocchio::JointIndex joint =
                    state_->get_pinocchio()->frames[data.frame].parent;
                pinocchio::getJointAccelerationDerivatives(
                    *state_->get_pinocchio().get(), *data.pinocchio, joint, pinocchio::LOCAL,
                    data.v_partial_dq, data.a_partial_dq, data.a_partial_dv, data.a_partial_da);
                const std::size_t nv = state_->get_nv();
                pinocchio::skew(data.v.linear(), data.vv_skew);
                pinocchio::skew(data.v.angular(), data.vw_skew);
                data.fXjdv_dq.noalias() = data.fXj * data.v_partial_dq;
                data.fXjda_dq.noalias() = data.fXj * data.a_partial_dq;
                data.fXjda_dv.noalias() = data.fXj * data.a_partial_dv;
                data.da0_local_dx.leftCols(nv) = data.fXjda_dq.template topRows<3>();
                data.da0_local_dx.leftCols(nv).noalias() +=
                    data.vw_skew * data.fXjdv_dq.template topRows<3>();
                data.da0_local_dx.leftCols(nv).noalias() -=
                    data.vv_skew * data.fXjdv_dq.template bottomRows<3>();
                data.da0_local_dx.rightCols(nv) = data.fXjda_dv.template topRows<3>();
                data.da0_local_dx.rightCols(nv).noalias() +=
                    data.vw_skew * data.fJf.template topRows<3>();
                data.da0_local_dx.rightCols(nv).noalias() -=
                    data.vv_skew * data.fJf.template bottomRows<3>();
                const Eigen::Ref<const Matrix3s> oRf = data.pinocchio->oMf[id_].rotation();

                if (gains_[0] != 0.)
                {
                    pinocchio::skew(data.dp_local, data.dp_skew);
                    data.da0_local_dx.leftCols(nv).noalias() +=
                        gains_[0] * data.dp_skew * data.fJf.template bottomRows<3>();
                    data.da0_local_dx.leftCols(nv).noalias() +=
                        gains_[0] * data.fJf.template topRows<3>();
                }
                if (gains_[1] != 0.)
                {
                    data.da0_local_dx.leftCols(nv).noalias() +=
                        gains_[1] * data.fXjdv_dq.template topRows<3>();
                    data.da0_local_dx.rightCols(nv).noalias() +=
                        gains_[1] * data.fJf.template topRows<3>();
                }
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    data.da0_dx.row(0) = (Raxis_ * data.da0_local_dx).row(2);
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    // Recalculate the constrained accelerations after imposing contact
                    // constraints. This is necessary for the forward-dynamics case.
                    data.a0_local = pinocchio::getFrameClassicalAcceleration(
                                        *state_->get_pinocchio().get(), *data.pinocchio, id_,
                                        pinocchio::LOCAL)
                                        .linear();
                    if (gains_[0] != 0.)
                    {
                        data.a0_local += gains_[0] * data.dp_local;
                    }
                    if (gains_[1] != 0.)
                    {
                        data.a0_local += gains_[1] * data.v.linear();
                    }
                    data.a0[0] = (Raxis_ * oRf * data.a0_local)[2];

                    pinocchio::skew((Raxis_ * oRf * data.a0_local).template head<3>(),
                                    data.a0_skew);
                    data.a0_world_skew.noalias() = data.a0_skew * Raxis_ * oRf;
                    data.da0_dx.row(0) = (Raxis_ * oRf * data.da0_local_dx).row(2);
                    data.da0_dx.leftCols(nv).row(0) -=
                        (data.a0_world_skew * data.fJf.template bottomRows<3>()).row(2);
                    break;
                }
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataDerived &data,
                             const Eigen::MatrixBase<ForceVectorType> &f)
            {
                const Eigen::Ref<const Matrix3s> R = data.jMf.rotation();
                data.f.linear()[2] = force[0];
                data.f.linear().template head<2>().setZero();
                data.f.angular().setZero();
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    data.fext.linear() = (R * Raxis_.transpose()).col(2) * force[0];
                    data.fext.angular() = data.jMf.translation().cross(data.fext.linear());
                    data.dtau_dq.setZero();
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    const Eigen::Ref<const Matrix3s> oRf = data.pinocchio->oMf[id_].rotation();
                    data.f_local.linear().noalias() =
                        (oRf.transpose() * Raxis_.transpose()).col(2) * force[0];
                    data.f_local.angular().setZero();
                    data.fext = data.jMf.act(data.f_local);
                    pinocchio::skew(data.f_local.linear(), data.f_skew);
                    data.fJf_df.noalias() = data.f_skew * data.fJf.template bottomRows<3>();
                    data.dtau_dq.noalias() =
                        -data.fJf.template topRows<3>().transpose() * data.fJf_df;
                    break;
                }
            }

        }; // struct ContactModel1dTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_1d_hpp__
