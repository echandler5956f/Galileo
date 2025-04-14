#ifndef __galileo_multibody_contacts_contact_2d_hpp__
#define __galileo_multibody_contacts_contact_2d_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{

    template <typename RobotSpec>
    struct Contact2dTpl;

    template <typename RobotSpec>
    struct traits<Contact2dTpl<RobotSpec>>
    {
        using RS = RobotSpec;

        using VarScalar = typename RS::VarScalar;
        using NumScalar = typename RS::NumScalar;
        static constexpr int Options = RS::Options;

        using ContactDataDerived = ContactData2dTpl<RobotSpec>;
        using ContactModelDerived = ContactModel2dTpl<RobotSpec>;

        static constexpr int NC = 1;
        // using RobotData_t = // TODO: add robot data
        using MatrixNcNv_t = Eigen::Matrix<VarScalar, NC, RS::NV, Options>;
        using MatrixNcNdx_t = Eigen::Matrix<VarScalar, NC, RS::NDX, Options>;
        using MatrixNcNu_t = Eigen::Matrix<VarScalar, NC, RS::NU, Options>;
    };

    template <typename RobotSpec>
    struct traits<ContactData2dTpl<RobotSpec>>
    {
        using ContactDerived = Contact2dTpl<RobotSpec>;
        using VarScalar = traits<ContactDerived>::VarScalar;
        using NumScalar = traits<ContactDerived>::NumScalar;
    };

    template <typename RobotSpec>
    struct traits<ContactModel2dTpl<RobotSpec>>
    {
        using ContactDerived = Contact2dTpl<RobotSpec>;
        using VarScalar = traits<ContactDerived>::VarScalar;
        using NumScalar = traits<ContactDerived>::NumScalar;
    };

    template <typename RobotSpec>
    struct ContactData2dTpl : ContactDataBase<ContactData2dTpl<RobotSpec>, RobotSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

    template <typename RobotSpec>
    struct ContactModel2dTpl : ContactDataBase<ContactModel2dTpl<RobotSpec>, RobotSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using RS = RobotSpec;

        using ContactDerived = Contact2dTpl<RobotSpec>;
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
            d->a = pinocchio::getFrameAcceleration(*state_->get_pinocchio().get(),
                                                   *d->pinocchio, id_);

            d->Jc.row(0) = d->fJf.row(0);
            d->Jc.row(1) = d->fJf.row(2);

            d->vw = d->v.angular();
            d->vv = d->v.linear();

            d->a0[0] = d->a.linear()[0] + d->vw[1] * d->vv[2] - d->vw[2] * d->vv[1];
            d->a0[1] = d->a.linear()[2] + d->vw[0] * d->vv[1] - d->vw[1] * d->vv[0];

            if (gains_[0] != 0.)
            {
                d->a0[0] +=
                    gains_[0] * (d->pinocchio->oMf[id_].translation()[0] - xref_[0]);
                d->a0[1] +=
                    gains_[0] * (d->pinocchio->oMf[id_].translation()[2] - xref_[1]);
            }
            if (gains_[1] != 0.)
            {
                d->a0[0] += gains_[1] * d->vv[0];
                d->a0[1] += gains_[1] * d->vv[2];
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
            pinocchio::skew(d->vv, d->vv_skew);
            pinocchio::skew(d->vw, d->vw_skew);
            d->fXjdv_dq.noalias() = d->fXj * d->v_partial_dq;
            d->fXjda_dq.noalias() = d->fXj * d->a_partial_dq;
            d->fXjda_dv.noalias() = d->fXj * d->a_partial_dv;

            d->da0_dx.leftCols(nv).row(0) = d->fXjda_dq.row(0);
            d->da0_dx.leftCols(nv).row(0).noalias() +=
                d->vw_skew.row(0) * d->fXjdv_dq.template topRows<3>();
            d->da0_dx.leftCols(nv).row(0).noalias() -=
                d->vv_skew.row(0) * d->fXjdv_dq.template bottomRows<3>();

            d->da0_dx.leftCols(nv).row(1) = d->fXjda_dq.row(2);
            d->da0_dx.leftCols(nv).row(1).noalias() +=
                d->vw_skew.row(2) * d->fXjdv_dq.template topRows<3>();
            d->da0_dx.leftCols(nv).row(1).noalias() -=
                d->vv_skew.row(2) * d->fXjdv_dq.template bottomRows<3>();

            d->da0_dx.rightCols(nv).row(0) = d->fXjda_dv.row(0);
            d->da0_dx.rightCols(nv).row(0).noalias() +=
                d->vw_skew.row(0) * d->fJf.template topRows<3>();
            d->da0_dx.rightCols(nv).row(0).noalias() -=
                d->vv_skew.row(0) * d->fJf.template bottomRows<3>();

            d->da0_dx.rightCols(nv).row(1) = d->fXjda_dv.row(2);
            d->da0_dx.rightCols(nv).row(1).noalias() +=
                d->vw_skew.row(2) * d->fJf.template topRows<3>();
            d->da0_dx.rightCols(nv).row(1).noalias() -=
                d->vv_skew.row(2) * d->fJf.template bottomRows<3>();

            if (gains_[0] != 0.)
            {
                const Eigen::Ref<const Matrix3s> oRf = d->pinocchio->oMf[id_].rotation();
                d->oRf(0, 0) = oRf(0, 0);
                d->oRf(1, 0) = oRf(2, 0);
                d->oRf(0, 1) = oRf(0, 2);
                d->oRf(1, 1) = oRf(2, 2);
                d->da0_dx.leftCols(nv).noalias() += gains_[0] * d->oRf * d->Jc;
            }
            if (gains_[1] != 0.)
            {
                d->da0_dx.leftCols(nv).row(0).noalias() +=
                    gains_[1] * d->fXj.row(0) * d->v_partial_dq;
                d->da0_dx.leftCols(nv).row(1).noalias() +=
                    gains_[1] * d->fXj.row(2) * d->v_partial_dq;
                d->da0_dx.rightCols(nv).row(0).noalias() +=
                    gains_[1] * d->fXj.row(0) * d->a_partial_da;
                d->da0_dx.rightCols(nv).row(1).noalias() +=
                    gains_[1] * d->fXj.row(2) * d->a_partial_da;
            }
        }

        template <typename ForceVectorType>
        void updateForce(ContactDataDerived &data,
                         const Eigen::MatrixBase<ForceVectorType> &f)
        {
            Data *d = static_cast<Data *>(data.get());
            const Eigen::Ref<const Matrix3s> R = d->jMf.rotation();
            data->f.linear() = R.col(0) * force[0] + R.col(2) * force[1];
            data->f.angular().setZero();
            data->fext.linear() = R.col(0) * force[0] + R.col(2) * force[1];
            data->fext.angular() = d->jMf.translation().cross(data->fext.linear());
        }

    }; // struct ContactModel2dTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_2d_hpp__
