#ifndef __galileo_multibody_contacts_contact_3d_hpp__
#define __galileo_multibody_contacts_contact_3d_hpp__

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/multibody/data.hpp>
#include <pinocchio/spatial/motion.hpp>

#include "galileo/multibody/contacts/contact-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct Contact3dTpl;

    template <typename PhaseSpec>
    struct traits<Contact3dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = ContactModel3dTpl<PS>;
        using Data_t = ContactData3dTpl<PS>;

        static constexpr int NC = 3;
        static constexpr int NU = traits<typename PS::NodeMeta_t>::NU;

        // Traits required by ForceDataBase
        using MatrixNcNv_t = Eigen::Matrix<VarScalar, NC, NV, Options>;
        using MatrixNcNdx_t = Eigen::Matrix<VarScalar, NC, NDX, Options>;
        using MatrixNcNu_t = Eigen::Matrix<VarScalar, NC, NU, Options>;

        // Traits required by ContactDataBase
        using VectorNc_t = Eigen::Matrix<VarScalar, NC, 1, Options>;
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
    struct ContactData3dTpl : public ContactDataBase<ContactData3dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        using Meta_t = Contact3dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

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

        // Accessor implementations required by ContactDataBase
        DEFAULT_ACCESSOR(ActionMatrix_t, fXj);
        DEFAULT_ACCESSOR(VectorNc_t, a0);
        DEFAULT_ACCESSOR(MatrixNcNdx_t, da0_dx);
        DEFAULT_ACCESSOR(MatrixNv_t, dtau_dq);

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

        // Members required by ContactDataBase
        ActionMatrix_t fXj;
        VectorNc_t a0;
        MatrixNcNdx_t da0_dx;
        MatrixNv_t dtau_dq;

        // Members used for ContactModel3dTpl
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

        template <typename DataCollector>
        ContactData3dTpl(const Model_t &model, DataCollector *const collector)
            : robot(collector->robot),
              frame(0),
              type(model.type()),
              jMf(SE3_t::Identity()),
              Jc(3, NV),
              f(Force_t::Zero()),
              fext(Force_t::Zero()),
              df_dx(3, NDX),
              df_du(3, model.nu()),
              fXj(jMf.inverse().toActionMatrix()),
              a0(3),
              da0_dx(3, NDX),
              dtau_dq(NV, NV),
              v(Motion_t::Zero()),
              f_local(Force_t::Zero()),
              da0_local_dx(3, NDX),
              fJf(6, NV),
              v_partial_dq(6, NV),
              a_partial_dq(6, NV),
              a_partial_dv(6, NV),
              a_partial_da(6, NV),
              fXjdv_dq(6, NV),
              fXjda_dq(6, NV),
              fXjda_dv(6, NV),
              fJf_df(3, NV)
        {
            Jc.setZero();
            df_dx.setZero();
            df_du.setZero();
            a0.setZero();
            da0_dx.setZero();
            dtau_dq.setZero();
            frame = model.id();
            jMf = model.robot()->frames[frame].placement;
            fXj = jMf.inverse().toActionMatrix();
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
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = Contact3dTpl<PhaseSpec>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Base_t = ContactModelBase<ContactModel3dTpl<PhaseSpec>, PhaseSpec>;
        using Base_t::updateForceDiff;
        using Base_t::setZeroForce;
        using Base_t::setZeroForceDiff;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        ContactModel3dTpl(State_t *state, const FrameIndex_t id,
                          const Vector3_t &xref, const ReferenceFrame_t type,
                          const int nu, const Vector2_t &gains) : robot_(state->get_robot()), id_(id), xref_(xref), type_(type), nu_(nu), gains_(gains)
        {
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            auto q = x.template head<PS::NQ>();
            auto v = x.template tail<PS::NV>();

            pinocchio::updateFramePlacement(*robot_, *(data.robot),
                                            id_);
            pinocchio::getFrameJacobian(*robot_, *(data.robot),
                                        id_, pinocchio::LOCAL, data.fJf);
            data.v = pinocchio::getFrameVelocity(*robot_,
                                                 *(data.robot), id_);
            data.a0_local =
                pinocchio::getFrameClassicalAcceleration(
                    *robot_, *(data.robot), id_, pinocchio::LOCAL)
                    .linear();

            const Eigen::Ref<const Matrix3_t> oRf = data.robot->oMf[id_].rotation();
            if (gains_[0] != 0.)
            {
                data.dp = data.robot->oMf[id_].translation() - xref_;
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
                data.Jc = data.fJf.template topRows<3>();
                data.a0 = data.a0_local;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                data.Jc.noalias() = oRf * data.fJf.template topRows<3>();
                data.a0.noalias() = oRf * data.a0_local;
                break;
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const pinocchio::JointIndex joint =
                robot_->frames[data.frame].parent;
            pinocchio::getJointAccelerationDerivatives(
                *robot_, *(data.robot), joint, pinocchio::LOCAL,
                data.v_partial_dq, data.a_partial_dq, data.a_partial_dv, data.a_partial_da);
            pinocchio::skew(data.v.linear(), data.vv_skew);
            pinocchio::skew(data.v.angular(), data.vw_skew);
            data.fXjdv_dq.noalias() = data.fXj * data.v_partial_dq;
            data.fXjda_dq.noalias() = data.fXj * data.a_partial_dq;
            data.fXjda_dv.noalias() = data.fXj * data.a_partial_dv;
            data.da0_local_dx.leftCols(NV) = data.fXjda_dq.template topRows<3>();
            data.da0_local_dx.leftCols(NV).noalias() +=
                data.vw_skew * data.fXjdv_dq.template topRows<3>();
            data.da0_local_dx.leftCols(NV).noalias() -=
                data.vv_skew * data.fXjdv_dq.template bottomRows<3>();
            data.da0_local_dx.rightCols(NV) = data.fXjda_dv.template topRows<3>();
            data.da0_local_dx.rightCols(NV).noalias() +=
                data.vw_skew * data.fJf.template topRows<3>();
            data.da0_local_dx.rightCols(NV).noalias() -=
                data.vv_skew * data.fJf.template bottomRows<3>();
            const Eigen::Ref<const Matrix3_t> oRf = data.robot->oMf[id_].rotation();

            if (gains_[0] != 0.)
            {
                pinocchio::skew(data.dp_local, data.dp_skew);
                data.da0_local_dx.leftCols(NV).noalias() +=
                    gains_[0] * data.dp_skew * data.fJf.template bottomRows<3>();
                data.da0_local_dx.leftCols(NV).noalias() +=
                    gains_[0] * data.fJf.template topRows<3>();
            }
            if (gains_[1] != 0.)
            {
                data.da0_local_dx.leftCols(NV).noalias() +=
                    gains_[1] * data.fXjdv_dq.template topRows<3>();
                data.da0_local_dx.rightCols(NV).noalias() +=
                    gains_[1] * data.fJf.template topRows<3>();
            }
            switch (type_)
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.da0_dx = data.da0_local_dx;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                // Recalculate the constrained accelerations after imposing contact
                // constraints. This is necessary for the forward-dynamics case.
                data.a0_local = pinocchio::getFrameClassicalAcceleration(
                                    *robot_, *(data.robot), id_,
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
                data.a0.noalias() = oRf * data.a0_local;

                pinocchio::skew(data.a0.template head<3>(), data.a0_skew);
                data.a0_world_skew.noalias() = data.a0_skew * oRf;
                data.da0_dx.noalias() = oRf * data.da0_local_dx;
                data.da0_dx.leftCols(NV).noalias() -=
                    data.a0_world_skew * data.fJf.template bottomRows<3>();
                break;
            }
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            data.f.linear() = force;
            data.f.angular().setZero();
            switch (type_)
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.fext = data.jMf.act(data.f);
                data.dtau_dq.setZero();
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                const Eigen::Ref<const Matrix3_t> oRf = data.robot->oMf[id_].rotation();
                data.f_local.linear().noalias() = oRf.transpose() * force;
                data.f_local.angular().setZero();
                data.fext = data.jMf.act(data.f_local);
                pinocchio::skew(data.f_local.linear(), data.f_skew);
                data.fJf_df.noalias() = data.f_skew * data.fJf.template bottomRows<3>();
                data.dtau_dq.noalias() =
                    -data.fJf.template topRows<3>().transpose() * data.fJf_df;
                break;
            }
        }

        const RobotModel_t *robot_impl() const
        {
            return robot_;
        }

        FrameIndex_t id_impl() const
        {
            return id_;
        }

        void set_id_impl(const FrameIndex_t &id)
        {
            id_ = id;
        }

        ReferenceFrame_t type_impl() const
        {
            return type_;
        }

        void set_type_impl(const ReferenceFrame_t &type)
        {
            type_ = type;
        }

        int nu_impl() const
        {
            return nu_;
        }

        int nc_impl() const
        {
            return traits<Meta_t>::NC;
        }

    protected:
        const RobotModel_t *robot_;
        FrameIndex_t id_;
        ReferenceFrame_t type_;
        Vector3_t xref_;
        Vector2_t gains_;

        int nu_;

    }; // struct ContactModel3dTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_3d_hpp__
