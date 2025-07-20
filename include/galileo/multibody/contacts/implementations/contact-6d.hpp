#ifndef __galileo_multibody_contacts_contact_6d_hpp__
#define __galileo_multibody_contacts_contact_6d_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct Contact6dTpl;

    template <typename PhaseSpec>
    struct traits<Contact6dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = Contact6dTpl<PS>;
        using Model_t = ContactModel6dTpl<PS>;
        using Data_t = ContactData6dTpl<PS>;

        using DimNC_t = DimensionTpl<6>;

        // Traits required by ForceDataBase
        using MatrixNcNv_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNV_t::Value, PS::Options>;
        using MatrixNcNdx_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNU_t::Value, PS::Options>;

        // Traits required by ContactDataBase
        using VectorNc_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, 1, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ContactData6dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = Contact6dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ContactModel6dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = Contact6dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ContactData6dTpl
        : ContactDataBase<ContactData6dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = Contact6dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ContactDataBase<ContactData6dTpl<PS>, PS>;

        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        // Members required by ForceDataBase
        RobotData_t *robot;
        FrameIndex_t frame;
        ReferenceFrame_t type;
        SE3_t jMf;
        MatrixNcNv_t Jc;
        Force_t f;
        Force_t fext;
        MatrixNcNdx_t df_dx;
        MatrixNcNu_t df_du;

        // Accessor implementations required by ForceDataBase
        DEFAULT_ACCESSOR(RobotData_t *, robot);
        DEFAULT_ACCESSOR(FrameIndex_t, frame);
        DEFAULT_ACCESSOR(ReferenceFrame_t, type);
        DEFAULT_ACCESSOR(SE3_t, jMf);
        DEFAULT_ACCESSOR(MatrixNcNv_t, Jc);
        DEFAULT_ACCESSOR(Force_t, f);
        DEFAULT_ACCESSOR(Force_t, fext);
        DEFAULT_ACCESSOR(MatrixNcNdx_t, df_dx);
        DEFAULT_ACCESSOR(MatrixNcNu_t, df_du);

        // Members required by ContactDataBase
        ActionMatrix_t fXj;
        VectorNc_t a0;
        MatrixNcNdx_t da0_dx;
        MatrixNv_t dtau_dq;

        // Accessor implementations required by ContactDataBase
        DEFAULT_ACCESSOR(ActionMatrix_t, fXj);
        DEFAULT_ACCESSOR(VectorNc_t, a0);
        DEFAULT_ACCESSOR(MatrixNcNdx_t, da0_dx);
        DEFAULT_ACCESSOR(MatrixNv_t, dtau_dq);

        // Members used for ContactModel6dTpl
        // Notice that we do not need to expose accessors for these because they are specific to the 6D contact model
        SE3_t rMf;
        SE3_t lwaMl;
        Motion_t v;
        Motion_t a0_local;
        Force_t f_local;
        Matrix6Ndx_t da0_local_dx;
        Matrix6Nv_t fJf;
        Matrix6Nv_t v_partial_dq;
        Matrix6Nv_t a_partial_dq;
        Matrix6Nv_t a_partial_dv;
        Matrix6Nv_t a_partial_da;
        Matrix3_t av_world_skew;
        Matrix3_t aw_world_skew;
        Matrix3_t av_skew;
        Matrix3_t aw_skew;
        Matrix3_t fv_skew;
        Matrix3_t fw_skew;
        Matrix6_t rMf_Jlog6;
        Matrix6Nv_t fJf_df;

        template <typename DataCollector>
        ContactData6dTpl(const Model_t &model, DataCollector *const collector)
            : robot(collector->robot),
              frame(0),
              type(model.get_type()),
              jMf(SE3_t::Identity()),
              Jc(model.get_nc(), model.get_ps().get_nv()),
              f(Force_t::Zero()),
              fext(Force_t::Zero()),
              df_dx(model.get_nc(), model.get_ps().get_ndx()),
              df_du(model.get_nc(), model.get_ps().get_nu()),
              fXj(jMf.inverse().toActionMatrix()),
              a0(model.get_nc()),
              da0_dx(model.get_nc(), model.get_ps().get_ndx()),
              dtau_dq(model.get_ps().get_nv(), model.get_ps().get_nv()),
              v(Motion_t::Zero()),
              f_local(Force_t::Zero()),
              da0_local_dx(model.get_nc(), model.get_ps().get_ndx()),
              fJf(6, model.get_ps().get_nv()),
              v_partial_dq(6, model.get_ps().get_nv()),
              a_partial_dq(6, model.get_ps().get_nv()),
              a_partial_dv(6, model.get_ps().get_nv()),
              a_partial_da(6, model.get_ps().get_nv()),
              fXjdv_dq(6, model.get_ps().get_nv()),
              fXjda_dq(6, model.get_ps().get_nv()),
              fXjda_dv(6, model.get_ps().get_nv()),
              fJf_df(model.get_nc(), model.get_ps().get_nv())
        {
            Jc.setZero();
            df_dx.setZero();
            df_du.setZero();
            a0.setZero();
            da0_dx.setZero();
            dtau_dq.setZero();
            frame = model.get_id();
            jMf = model.get_robot().frames[frame].placement;
            fXj = jMf.inverse().toActionMatrix();
            da0_local_dx.setZero();
            fJf.setZero();
            v_partial_dq.setZero();
            a_partial_dq.setZero();
            a_partial_dv.setZero();
            a_partial_da.setZero();
            av_world_skew.setZero();
            aw_world_skew.setZero();
            av_skew.setZero();
            aw_skew.setZero();
            fv_skew.setZero();
            fw_skew.setZero();
            rMf_Jlog6.setZero();
            fJf_df.setZero();
        }
    };

    template <typename PhaseSpec>
    struct ContactModel6dTpl
        : ContactModelBase<ContactModel6dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = Contact6dTpl<PhaseSpec>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ContactModelBase<ContactModel6dTpl<PS>, PS>;

        using DimNC_t = typename traits<Meta_t>::DimNC_t;

        ContactModel6dTpl(const PS &ps,
                          const RobotModel_t &robot,
                          const FrameIndex_t id,
                          const ReferenceFrame_t &type,
                          const SE3_t &pref,
                          const Vector2_t &gains)
            : Base(ps, id, type, DimNC_t()),
              robot_(robot),
              pref_(pref),
              gains_(gains)
        {
        }

        template <typename StateVectorType>
        void calc(ContactDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            pinocchio::updateFramePlacement(get_robot(),
                                            *data.robot, get_id());
            pinocchio::getFrameJacobian(get_robot(), *data.robot,
                                        get_id(), pinocchio::LOCAL, data.fJf);
            data.a0_local = pinocchio::getFrameAcceleration(get_robot(),
                                                            *data.robot, get_id());

            if (gains_[0] != 0.)
            {
                data.rMf = pref_.actInv(data.robot->oMf[get_id()]);
                data.a0_local += gains_[0] * pinocchio::log6(data.rMf);
            }
            if (gains_[1] != 0.)
            {
                data.v = pinocchio::getFrameVelocity(get_robot(),
                                                     *data.robot, get_id());
                data.a0_local += gains_[1] * data.v;
            }
            switch (get_type())
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.Jc = data.fJf;
                data.a0 = data.a0_local.toVector();
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                data.lwaMl.rotation(data.robot->oMf[get_id()].rotation());
                data.Jc.noalias() = data.lwaMl.toActionMatrix() * data.fJf;
                data.a0.noalias() = data.lwaMl.act(data.a0_local).toVector();
                break;
            }
        }

        template <typename StateVectorType>
        void calcDiff(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const pinocchio::JointIndex joint =
                get_robot().frames[data.frame].parent;
            pinocchio::getJointAccelerationDerivatives(
                get_robot(), *data.robot, joint, pinocchio::LOCAL,
                data.v_partial_dq, data.a_partial_dq, data.a_partial_dv, data.a_partial_da);
            leftCols(data.da0_local_dx, get_ps().get_nv_dim()).noalias() = data.fXj * data.a_partial_dq;
            rightCols(data.da0_local_dx, get_ps().get_nv_dim()).noalias() = data.fXj * data.a_partial_dv;

            if (gains_[0] != 0.)
            {
                pinocchio::Jlog6(data.rMf, data.rMf_Jlog6);
                leftCols(data.da0_local_dx, get_ps().get_nv_dim()).noalias() += gains_[0] * data.rMf_Jlog6 * data.fJf;
            }
            if (gains_[1] != 0.)
            {
                leftCols(data.da0_local_dx, get_ps().get_nv_dim()).noalias() +=
                    gains_[1] * data.fXj * data.v_partial_dq;
                rightCols(data.da0_local_dx, get_ps().get_nv_dim()).noalias() += gains_[1] * data.fJf;
            }
            switch (get_type())
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.da0_dx = data.da0_local_dx;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                // Recalculate the constrained accelerations after imposing contact
                // constraints. This is necessary for the forward-dynamics case.
                data.a0_local = pinocchio::getFrameAcceleration(
                    get_robot(), *data.robot, get_id());
                if (gains_[0] != 0.)
                {
                    data.a0_local += gains_[0] * pinocchio::log6(data.rMf);
                }
                if (gains_[1] != 0.)
                {
                    data.a0_local += gains_[1] * data.v;
                }
                data.a0.noalias() = data.lwaMl.act(data.a0_local).toVector();

                const Eigen::Ref<const Matrix3_t> oRf = data.robot->oMf[get_id()].rotation();
                pinocchio::skew(head(data.a0, 3), data.av_skew);
                pinocchio::skew(tail(data.a0, 3), data.aw_skew);
                data.av_world_skew.noalias() = data.av_skew * oRf;
                data.aw_world_skew.noalias() = data.aw_skew * oRf;
                data.da0_dx.noalias() = data.lwaMl.toActionMatrix() * data.da0_local_dx;
                topRows(leftCols(data.da0_dx, get_ps().get_nv_dim()), 3).noalias() -=
                    data.av_world_skew * bottomRows(data.fJf, 3);
                bottomRows(leftCols(data.da0_dx, get_ps().get_nv_dim()), 3).noalias() -=
                    data.aw_world_skew * bottomRows(data.fJf, 3);
                break;
            }
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        template <typename ForceVectorType>
        void updateForce(ContactDataDerived &data,
                         const Eigen::MatrixBase<ForceVectorType> &f) const
        {
            data.f = pinocchio::ForceTpl<typename PS::VarScalar>(f);
            switch (get_type())
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.fext = data.jMf.act(data.f);
                data.dtau_dq.setZero();
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                data.f_local = data.lwaMl.actInv(data.f);
                data.fext = data.jMf.act(data.f_local);
                pinocchio::skew(data.f_local.linear(), data.fv_skew);
                pinocchio::skew(data.f_local.angular(), data.fw_skew);
                topRows(data.fJf_df, 3).noalias() =
                    data.fv_skew * bottomRows(data.fJf, 3);
                bottomRows(data.fJf_df, 3).noalias() =
                    data.fw_skew * bottomRows(data.fJf, 3);
                data.dtau_dq.noalias() = -data.fJf.transpose() * data.fJf_df;
                break;
            }
        }
        using Base::setZeroForce;
        using Base::setZeroForceDiff;
        using Base::updateForceDiff;

        const RobotModel_t &get_robot() const
        {
            return robot_.get();
        }

        using Base::get_ps;

        using Base::get_id;
        using Base::get_type;

        using Base::set_id;
        using Base::set_type;

        using Base::get_nc;
        using Base::get_nc_dim;

    protected:
        std::reference_wrapper<const RobotModel_t> robot_;
        SE3_t pref_;
        Vector2_t gains_;

    }; // struct ContactModel6dTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_6d_hpp__
