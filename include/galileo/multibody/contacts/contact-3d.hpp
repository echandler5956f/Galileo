#ifndef __galileo_multibody_contacts_contact_3d_hpp__
#define __galileo_multibody_contacts_contact_3d_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct Contact3dTpl;

    template <typename PhaseSpec>
    struct traits<Contact3dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_PINOCCIO_TYPES_TYPEDEF(PS);

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = ContactModel3dTpl<PS>;
        using Data_t = ContactData3dTpl<PS>;

        static constexpr int NC = 3;
        using MatrixNcNv_t = Eigen::Matrix<VarScalar, NC, PS::NV, Options>;
        using MatrixNcNdx_t = Eigen::Matrix<VarScalar, NC, PS::NDX, Options>;
        using MatrixNcNu_t = Eigen::Matrix<VarScalar, NC, PS::NU, Options>;
    };

    template <typename PhaseSpec>
    struct traits<ContactData3dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ContactModel3dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ContactData3dTpl : ContactDataBase<ContactData3dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_PINOCCIO_TYPES_TYPEDEF(PS);

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(RobotData_t *, robot_data_pointer);
        DEFAULT_ACCESSOR(Index_t, frame);
        DEFAULT_ACCESSOR(ReferenceFrame_t, type);
        DEFAULT_ACCESSOR(SE3_t, jMf);
        DEFAULT_ACCESSOR(MatrixNcNv_t, Jc);
        DEFAULT_ACCESSOR(Force_t, f);
        DEFAULT_ACCESSOR(Force_t, fext);
        DEFAULT_ACCESSOR(MatrixNcNdx_t, df_dx);
        DEFAULT_ACCESSOR(MatrixNcNu_t, df_du);

        DEFAULT_ACCESSOR(ActionMatrix_t, fXj);
        DEFAULT_ACCESSOR(VectorNc_t, a0);
        DEFAULT_ACCESSOR(MatrixNcNdx_t, da0_dx);
        DEFAULT_ACCESSOR(MatrixNv_t, dtau_dq);

        RobotData_t *robot_data_pointer;
        Index_t frame;
        ReferenceFrame_t type;
        SE3_t jMf;
        MatrixNcNv_t Jc;
        Force_t f;
        Force_t fext;
        MatrixNcNdx_t df_dx;
        MatrixNcNu_t df_du;
        
        ActionMatrix_t fXj;
        VectorNc_t a0;
        MatrixNcNdx_t da0_dx;
        MatrixNv_t dtau_dq;

    };

    template <typename PhaseSpec>
    struct ContactModel3dTpl : ContactDataBase<ContactModel3dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using ContactDerived = Contact3dTpl<PhaseSpec>;
        using ContactModelDerived = traits<ContactDerived>::ContactModelDerived;
        using ContactDataDerived = traits<ContactDerived>::ContactDataDerived;

        template <typename StateVectorType>
        void calc(ContactDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x)
        {
            Data *d = static_cast<Data *>(data.get());
            pinocchio::updateFramePlacement(*state_->get_pinocchio().get(), *d->pinocchio,
                                            id_);
            pinocchio::getFrameJacobian(*state_->get_pinocchio().get(), *d->pinocchio,
                                        id_, pinocchio::LOCAL, d->fJf);
            d->v = pinocchio::getFrameVelocity(*state_->get_pinocchio().get(),
                                               *d->pinocchio, id_);
            d->a0_local =
                pinocchio::getFrameClassicalAcceleration(
                    *state_->get_pinocchio().get(), *d->pinocchio, id_, pinocchio::LOCAL)
                    .linear();

            const Eigen::Ref<const Matrix3s> oRf = d->pinocchio->oMf[id_].rotation();
            if (gains_[0] != 0.)
            {
                d->dp = d->pinocchio->oMf[id_].translation() - xref_;
                d->dp_local.noalias() = oRf.transpose() * d->dp;
                d->a0_local += gains_[0] * d->dp_local;
            }
            if (gains_[1] != 0.)
            {
                d->a0_local += gains_[1] * d->v.linear();
            }
            switch (type_)
            {
            case pinocchio::ReferenceFrame::LOCAL:
                d->Jc = d->fJf.template topRows<3>();
                d->a0 = d->a0_local;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                d->Jc.noalias() = oRf * d->fJf.template topRows<3>();
                d->a0.noalias() = oRf * d->a0_local;
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
            pinocchio::skew(d->v.linear(), d->vv_skew);
            pinocchio::skew(d->v.angular(), d->vw_skew);
            d->fXjdv_dq.noalias() = d->fXj * d->v_partial_dq;
            d->fXjda_dq.noalias() = d->fXj * d->a_partial_dq;
            d->fXjda_dv.noalias() = d->fXj * d->a_partial_dv;
            d->da0_local_dx.leftCols(nv) = d->fXjda_dq.template topRows<3>();
            d->da0_local_dx.leftCols(nv).noalias() +=
                d->vw_skew * d->fXjdv_dq.template topRows<3>();
            d->da0_local_dx.leftCols(nv).noalias() -=
                d->vv_skew * d->fXjdv_dq.template bottomRows<3>();
            d->da0_local_dx.rightCols(nv) = d->fXjda_dv.template topRows<3>();
            d->da0_local_dx.rightCols(nv).noalias() +=
                d->vw_skew * d->fJf.template topRows<3>();
            d->da0_local_dx.rightCols(nv).noalias() -=
                d->vv_skew * d->fJf.template bottomRows<3>();
            const Eigen::Ref<const Matrix3s> oRf = d->pinocchio->oMf[id_].rotation();

            if (gains_[0] != 0.)
            {
                pinocchio::skew(d->dp_local, d->dp_skew);
                d->da0_local_dx.leftCols(nv).noalias() +=
                    gains_[0] * d->dp_skew * d->fJf.template bottomRows<3>();
                d->da0_local_dx.leftCols(nv).noalias() +=
                    gains_[0] * d->fJf.template topRows<3>();
            }
            if (gains_[1] != 0.)
            {
                d->da0_local_dx.leftCols(nv).noalias() +=
                    gains_[1] * d->fXjdv_dq.template topRows<3>();
                d->da0_local_dx.rightCols(nv).noalias() +=
                    gains_[1] * d->fJf.template topRows<3>();
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
                d->a0_local = pinocchio::getFrameClassicalAcceleration(
                                  *state_->get_pinocchio().get(), *d->pinocchio, id_,
                                  pinocchio::LOCAL)
                                  .linear();
                if (gains_[0] != 0.)
                {
                    d->a0_local += gains_[0] * d->dp_local;
                }
                if (gains_[1] != 0.)
                {
                    d->a0_local += gains_[1] * d->v.linear();
                }
                d->a0.noalias() = oRf * d->a0_local;

                pinocchio::skew(d->a0.template head<3>(), d->a0_skew);
                d->a0_world_skew.noalias() = d->a0_skew * oRf;
                d->da0_dx.noalias() = oRf * d->da0_local_dx;
                d->da0_dx.leftCols(nv).noalias() -=
                    d->a0_world_skew * d->fJf.template bottomRows<3>();
                break;
            }
        }

        template <typename ForceVectorType>
        void updateForce(ContactDataDerived &data,
                         const Eigen::MatrixBase<ForceVectorType> &f)
        {
            Data *d = static_cast<Data *>(data.get());
            data->f.linear() = force;
            data->f.angular().setZero();
            switch (type_)
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data->fext = d->jMf.act(data->f);
                data->dtau_dq.setZero();
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                const Eigen::Ref<const Matrix3s> oRf = d->pinocchio->oMf[id_].rotation();
                d->f_local.linear().noalias() = oRf.transpose() * force;
                d->f_local.angular().setZero();
                data->fext = data->jMf.act(d->f_local);
                pinocchio::skew(d->f_local.linear(), d->f_skew);
                d->fJf_df.noalias() = d->f_skew * d->fJf.template bottomRows<3>();
                data->dtau_dq.noalias() =
                    -d->fJf.template topRows<3>().transpose() * d->fJf_df;
                break;
            }
        }

    }; // struct ContactModel3dTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_3d_hpp__
