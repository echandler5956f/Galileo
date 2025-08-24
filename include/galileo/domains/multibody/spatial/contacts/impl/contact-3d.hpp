#ifndef __galileo_multibody_spatial_contacts_contact_3d_hpp__
#define __galileo_multibody_spatial_contacts_contact_3d_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/multibody/data.hpp>

#include "galileo/domains/multibody/spatial/contacts/contact-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct Contact3dTpl;

    template <typename PhaseSpec>
    struct traits<Contact3dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = ContactModel3dTpl<PS>;
        using Data_t = ContactData3dTpl<PS>;

        using DimNC_t = DimensionTpl<3>;
        static constexpr int NC = DimNC_t::Value;

        // Traits required by ForceDataBase
        using MatrixNcNv_t = Eigen::GMatrix<typename PS::VarScalar, NC, PS::NV, PS::Options>;
        using MatrixNcNdx_t = Eigen::GMatrix<typename PS::VarScalar, NC, PS::NDX, PS::Options>;
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, NC, PS::NU, PS::Options>;

        // Traits required by ContactDataBase
        using VectorNc_t = Eigen::GMatrix<typename PS::VarScalar, NC, 1, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ContactData3dTpl<PhaseSpec>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = Contact3dTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ContactModel3dTpl<PhaseSpec>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = Contact3dTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ContactData3dTpl : public ContactDataBase<ContactData3dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ContactDataBase<ContactData3dTpl<PS>, PS>;

        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        using RobotData_t = typename PS::RobotData_t;
        using FrameIndex_t = typename PS::FrameIndex_t;
        using ReferenceFrame_t = typename PS::ReferenceFrame_t;
        using SE3_t = typename PS::SE3_t;
        using ActionMatrix_t = typename PS::ActionMatrix_t;
        using Force_t = typename PS::Force_t;
        using Motion_t = typename PS::Motion_t;
        using MatrixNv_t = typename PS::MatrixNv_t;
        using Vector3_t = typename PS::Vector3_t;
        using Matrix3_t = typename PS::Matrix3_t;
        using Matrix6X_t = typename PS::Matrix6X_t;
        using Matrix3X_t = typename PS::Matrix3X_t;

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

        // Members used for ContactModel3dTpl
        // Notice that we do not need to expose accessors for these because they are specific to the 3D contact model
        Motion_t v;
        Vector3_t a0_local;
        Vector3_t dp;
        Vector3_t dp_local;
        Force_t f_local;
        MatrixNcNdx_t da0_local_dx;
        Matrix6X_t fJf;
        Matrix6X_t v_partial_dq;
        Matrix6X_t a_partial_dq;
        Matrix6X_t a_partial_dv;
        Matrix6X_t a_partial_da;
        Matrix3_t vv_skew;
        Matrix3_t vw_skew;
        Matrix3_t a0_skew;
        Matrix3_t a0_world_skew;
        Matrix3_t dp_skew;
        Matrix3_t f_skew;
        Matrix6X_t fXjdv_dq;
        Matrix6X_t fXjda_dq;
        Matrix6X_t fXjda_dv;
        Matrix3X_t fJf_df;

        ContactData3dTpl(const Model_t &model, RobotData_t *const robot_data)
            : robot(robot_data),
              frame(model.get_id()),
              type(model.get_type()),
              jMf(model.get_state().get_robot().frames[frame].placement),
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
              a0_local(model.get_nc()),
              dp(model.get_nc()),
              dp_local(model.get_nc()),
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
            a0_local.setZero();
            dp.setZero();
            dp_local.setZero();
            da0_local_dx.setZero();
            fJf.setZero();
            v_partial_dq.setZero();
            a_partial_dq.setZero();
            a_partial_dv.setZero();
            a_partial_da.setZero();
            vv_skew.setZero();
            vw_skew.setZero();
            a0_skew.setZero();
            a0_world_skew.setZero();
            dp_skew.setZero();
            f_skew.setZero();
            fXjdv_dq.setZero();
            fXjda_dq.setZero();
            fXjda_dv.setZero();
            fJf_df.setZero();
        }
    };

    template <typename PhaseSpec>
    struct ContactModel3dTpl : public ContactModelBase<ContactModel3dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = Contact3dTpl<PhaseSpec>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ContactModelBase<ContactModel3dTpl<PS>, PS>;

        using DimNC_t = typename traits<Meta_t>::DimNC_t;

        using RobotData_t = typename PS::RobotData_t;
        using State_t = typename PS::State_t;
        using FrameIndex_t = typename PS::FrameIndex_t;
        using ReferenceFrame_t = typename PS::ReferenceFrame_t;
        using JointIndex_t = typename PS::JointIndex_t;
        using Vector3_t = typename PS::Vector3_t;
        using Vector2_t = typename PS::Vector2_t;

        template <typename Vector3Type, typename Vector2Type>
        ContactModel3dTpl(const PS &ps,
                          const State_t &state,
                          const FrameIndex_t id,
                          const ReferenceFrame_t &type,
                          const Eigen::MatrixBase<Vector3Type> &xref,
                          const Eigen::MatrixBase<Vector2Type> &gains)
            : Base(id, type, DimNC_t()), ps_(ps), state_(state), xref_(xref), gains_(gains)
        {
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            pinocchio::updateFramePlacement(get_state().get_robot(), *data.robot, get_id());
            pinocchio::getFrameJacobian(get_state().get_robot(), *data.robot, get_id(), pinocchio::LOCAL, data.fJf);
            data.v = pinocchio::getFrameVelocity(get_state().get_robot(), *data.robot, get_id());
            data.a0_local = pinocchio::getFrameClassicalAcceleration(
                                get_state().get_robot(), *data.robot, get_id(), pinocchio::LOCAL)
                                .linear();

            const auto oRf = data.robot->oMf[get_id()].rotation();
            if (gains_[0] != 0.)
            {
                data.dp = data.robot->oMf[get_id()].translation() - xref_;
                data.dp_local.noalias() = oRf.transpose() * data.dp;
                data.a0_local += gains_[0] * data.dp_local;
            }
            if (gains_[1] != 0.)
            {
                data.a0_local += gains_[1] * data.v.linear();
            }
            switch (get_type())
            {
            case ReferenceFrame_t::LOCAL:
                data.Jc = topRows<3>(data.fJf);
                data.a0 = data.a0_local;
                break;
            case ReferenceFrame_t::WORLD:
            case ReferenceFrame_t::LOCAL_WORLD_ALIGNED:
                data.Jc.noalias() = oRf * topRows<3>(data.fJf);
                data.a0.noalias() = oRf * data.a0_local;
                break;
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const auto nv_dim = get_ps().get_nv_dim();
            const JointIndex_t joint = get_state().get_robot().frames[data.frame].parentJoint;
            pinocchio::getJointAccelerationDerivatives(get_state().get_robot(),
                                                       *data.robot,
                                                       joint,
                                                       pinocchio::LOCAL,
                                                       data.v_partial_dq,
                                                       data.a_partial_dq,
                                                       data.a_partial_dv,
                                                       data.a_partial_da);
            pinocchio::skew(data.v.linear(), data.vv_skew);
            pinocchio::skew(data.v.angular(), data.vw_skew);
            data.fXjdv_dq.noalias() = data.fXj * data.v_partial_dq;
            data.fXjda_dq.noalias() = data.fXj * data.a_partial_dq;
            data.fXjda_dv.noalias() = data.fXj * data.a_partial_dv;
            leftCols(data.da0_local_dx, nv_dim) = topRows<3>(data.fXjda_dq);
            leftCols(data.da0_local_dx, nv_dim).noalias() += data.vw_skew * topRows<3>(data.fXjdv_dq);
            leftCols(data.da0_local_dx, nv_dim).noalias() -= data.vv_skew * bottomRows<3>(data.fXjdv_dq);
            rightCols(data.da0_local_dx, nv_dim) = topRows<3>(data.fXjda_dv);
            rightCols(data.da0_local_dx, nv_dim).noalias() += data.vw_skew * topRows<3>(data.fJf);
            rightCols(data.da0_local_dx, nv_dim).noalias() -= data.vv_skew * bottomRows<3>(data.fJf);
            const auto oRf = data.robot->oMf[get_id()].rotation();

            if (gains_[0] != 0.)
            {
                pinocchio::skew(data.dp_local, data.dp_skew);
                leftCols(data.da0_local_dx, nv_dim).noalias() += gains_[0] * data.dp_skew * bottomRows<3>(data.fJf);
                leftCols(data.da0_local_dx, nv_dim).noalias() += gains_[0] * topRows<3>(data.fJf);
            }
            if (gains_[1] != 0.)
            {
                leftCols(data.da0_local_dx, nv_dim).noalias() += gains_[1] * topRows<3>(data.fXjdv_dq);
                rightCols(data.da0_local_dx, nv_dim).noalias() += gains_[1] * topRows<3>(data.fJf);
            }
            switch (get_type())
            {
            case ReferenceFrame_t::LOCAL:
                data.da0_dx = data.da0_local_dx;
                break;
            case ReferenceFrame_t::WORLD:
            case ReferenceFrame_t::LOCAL_WORLD_ALIGNED:
                // Recalculate the constrained accelerations after imposing contact
                // constraints. This is necessary for the forward-dynamics case.
                data.a0_local = pinocchio::getFrameClassicalAcceleration(
                                    get_state().get_robot(), *data.robot, get_id(), pinocchio::LOCAL)
                                    .linear();
                if (gains_[0] != 0.)
                {
                    data.a0_local += gains_[0] * data.dp_local;
                }
                if (gains_[1] != 0.)
                {
                    data.a0_local += gains_[1] * data.v.linear();
                }
                data.a0.noalias() = oRf * data.a0_local;

                pinocchio::skew(head<3>(data.a0), data.a0_skew);
                data.a0_world_skew.noalias() = data.a0_skew * oRf;
                data.da0_dx.noalias() = oRf * data.da0_local_dx;
                leftCols(data.da0_dx, nv_dim).noalias() -= data.a0_world_skew * bottomRows<3>(data.fJf);
                break;
            }
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data, const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            data.f.linear() = force;
            data.f.angular().setZero();
            switch (get_type())
            {
            case ReferenceFrame_t::LOCAL:
                data.fext = data.jMf.act(data.f);
                data.dtau_dq.setZero();
                break;
            case ReferenceFrame_t::WORLD:
            case ReferenceFrame_t::LOCAL_WORLD_ALIGNED:
                const auto oRf = data.robot->oMf[get_id()].rotation();
                data.f_local.linear().noalias() = oRf.transpose() * force;
                data.f_local.angular().setZero();
                data.fext = data.jMf.act(data.f_local);
                pinocchio::skew(data.f_local.linear(), data.f_skew);
                data.fJf_df.noalias() = data.f_skew * bottomRows<3>(data.fJf);
                data.dtau_dq.noalias() = -topRows<3>(data.fJf).transpose() * data.fJf_df;
                break;
            }
        }

        Data_t createData(RobotData_t *const robot) const { return Data_t(*this, robot); }

        using Base::setZeroForce;
        using Base::setZeroForceDiff;
        using Base::updateForceDiff;
        using Base::get_id;
        using Base::get_type;
        using Base::set_id;
        using Base::set_type;
        using Base::get_nc;

        const PS &get_ps() const { return ps_; }
        const State_t &get_state() const { return state_; }

    protected:
        PS ps_;
        State_t state_;
        Vector3_t xref_;
        Vector2_t gains_;

    }; // struct ContactModel3dTpl

} // namespace galileo

#endif // __galileo_multibody_spatial_contacts_contact_3d_hpp__
