#ifndef __galileo_multibody_contacts_contact_6d_hpp__
#define __galileo_multibody_contacts_contact_6d_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename BasicSpec>
        struct Contact6dTpl;

        template <typename BasicSpec>
        struct traits<Contact6dTpl<BasicSpec>>
        {
            using BS = BasicSpec;

            using VarScalar = typename BS::VarScalar;
            using NumScalar = typename BS::NumScalar;
            static constexpr int Options = BS::Options;

            using ContactDataDerived = ContactData6dTpl<BasicSpec>;
            using ContactModelDerived = ContactModel6dTpl<BasicSpec>;

            static constexpr int NC = 1;
            // using RobotData_t = // TODO: add robot data
            using MatrixNcNv_t = Eigen::Matrix<VarScalar, NC, BS::NV, Options>;
            using MatrixNcNdx_t = Eigen::Matrix<VarScalar, NC, BS::NDX, Options>;
            using MatrixNcNu_t = Eigen::Matrix<VarScalar, NC, BS::NU, Options>;
        };

        template <typename BasicSpec>
        struct traits<ContactData6dTpl<BasicSpec>>
        {
            using ContactDerived = Contact6dTpl<BasicSpec>;
            using VarScalar = traits<ContactDerived>::VarScalar;
            using NumScalar = traits<ContactDerived>::NumScalar;
        };

        template <typename BasicSpec>
        struct traits<ContactModel6dTpl<BasicSpec>>
        {
            using ContactDerived = Contact6dTpl<BasicSpec>;
            using VarScalar = traits<ContactDerived>::VarScalar;
            using NumScalar = traits<ContactDerived>::NumScalar;
        };

        template <typename BasicSpec>
        struct ContactData6dTpl : ContactDataBase<ContactData6dTpl<BasicSpec>, BasicSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        template <typename BasicSpec>
        struct ContactModel6dTpl : ContactDataBase<ContactModel6dTpl<BasicSpec>, BasicSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using BS = BasicSpec;

            using ContactDerived = Contact6dTpl<BasicSpec>;
            using ContactModelDerived = traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = traits<ContactDerived>::ContactDataDerived;

            template <typename StateVectorType>
            void calc(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                Data *d = static_cast<Data *>(data.get());
                pinocchio::updateFramePlacement<Scalar>(*state_->get_pinocchio().get(),
                                                        *d->pinocchio, id_);
                pinocchio::getFrameJacobian(*state_->get_pinocchio().get(), *d->pinocchio,
                                            id_, pinocchio::LOCAL, d->fJf);
                d->a0_local = pinocchio::getFrameAcceleration(*state_->get_pinocchio().get(),
                                                              *d->pinocchio, id_);

                if (gains_[0] != 0.)
                {
                    d->rMf = pref_.actInv(d->pinocchio->oMf[id_]);
                    d->a0_local += gains_[0] * pinocchio::log6(d->rMf);
                }
                if (gains_[1] != 0.)
                {
                    d->v = pinocchio::getFrameVelocity(*state_->get_pinocchio().get(),
                                                       *d->pinocchio, id_);
                    d->a0_local += gains_[1] * d->v;
                }
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    data->Jc = d->fJf;
                    data->a0 = d->a0_local.toVector();
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    d->lwaMl.rotation(d->pinocchio->oMf[id_].rotation());
                    data->Jc.noalias() = d->lwaMl.toActionMatrix() * d->fJf;
                    data->a0.noalias() = d->lwaMl.act(d->a0_local).toVector();
                    break;
                }
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                Data *d = static_cast<Data *>(data.get());
                const pinocchio::JointIndex joint =
                    state_->get_pinocchio()->frames[d->frame].parent;
                pinocchio::getJointAccelerationDerivatives(
                    *state_->get_pinocchio().get(), *d->pinocchio, joint, pinocchio::LOCAL,
                    d->v_partial_dq, d->a_partial_dq, d->a_partial_dv, d->a_partial_da);
                const std::size_t nv = state_->get_nv();
                d->da0_local_dx.leftCols(nv).noalias() = d->fXj * d->a_partial_dq;
                d->da0_local_dx.rightCols(nv).noalias() = d->fXj * d->a_partial_dv;

                if (gains_[0] != 0.)
                {
                    pinocchio::Jlog6(d->rMf, d->rMf_Jlog6);
                    d->da0_local_dx.leftCols(nv).noalias() += gains_[0] * d->rMf_Jlog6 * d->fJf;
                }
                if (gains_[1] != 0.)
                {
                    d->da0_local_dx.leftCols(nv).noalias() +=
                        gains_[1] * d->fXj * d->v_partial_dq;
                    d->da0_local_dx.rightCols(nv).noalias() += gains_[1] * d->fJf;
                }
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    d->da0_dx = d->da0_local_dx;
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    // Recalculate the constrained accelerations after imposing contact
                    // constraints. This is necessary for the forward-dynamics case.
                    d->a0_local = pinocchio::getFrameAcceleration(
                        *state_->get_pinocchio().get(), *d->pinocchio, id_);
                    if (gains_[0] != 0.)
                    {
                        d->a0_local += gains_[0] * pinocchio::log6(d->rMf);
                    }
                    if (gains_[1] != 0.)
                    {
                        d->a0_local += gains_[1] * d->v;
                    }
                    data->a0.noalias() = d->lwaMl.act(d->a0_local).toVector();

                    const Eigen::Ref<const Matrix3s> oRf = d->pinocchio->oMf[id_].rotation();
                    pinocchio::skew(d->a0.template head<3>(), d->av_skew);
                    pinocchio::skew(d->a0.template tail<3>(), d->aw_skew);
                    d->av_world_skew.noalias() = d->av_skew * oRf;
                    d->aw_world_skew.noalias() = d->aw_skew * oRf;
                    d->da0_dx.noalias() = d->lwaMl.toActionMatrix() * d->da0_local_dx;
                    d->da0_dx.leftCols(nv).template topRows<3>().noalias() -=
                        d->av_world_skew * d->fJf.template bottomRows<3>();
                    d->da0_dx.leftCols(nv).template bottomRows<3>().noalias() -=
                        d->aw_world_skew * d->fJf.template bottomRows<3>();
                    break;
                }
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataDerived &data,
                             const Eigen::MatrixBase<ForceVectorType> &f)
            {
                if (force.size() != 6)
                {
                    throw_pretty("Invalid argument: "
                                 << "lambda has wrong dimension (it should be 6)");
                }
                Data *d = static_cast<Data *>(data.get());
                data->f = pinocchio::ForceTpl<Scalar>(force);
                switch (type_)
                {
                case pinocchio::ReferenceFrame::LOCAL:
                    data->fext = data->jMf.act(data->f);
                    data->dtau_dq.setZero();
                    break;
                case pinocchio::ReferenceFrame::WORLD:
                case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                    d->f_local = d->lwaMl.actInv(data->f);
                    data->fext = data->jMf.act(d->f_local);
                    pinocchio::skew(d->f_local.linear(), d->fv_skew);
                    pinocchio::skew(d->f_local.angular(), d->fw_skew);
                    d->fJf_df.template topRows<3>().noalias() =
                        d->fv_skew * d->fJf.template bottomRows<3>();
                    d->fJf_df.template bottomRows<3>().noalias() =
                        d->fw_skew * d->fJf.template bottomRows<3>();
                    d->dtau_dq.noalias() = -d->fJf.transpose() * d->fJf_df;
                    break;
                }
            }

        }; // struct ContactModel6dTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_6d_hpp__
