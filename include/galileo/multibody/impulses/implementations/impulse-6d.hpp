#ifndef __galileo_multibody_impulses_impulse_6d_hpp__
#define __galileo_multibody_impulses_impulse_6d_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/multibody/data.hpp>

#include "galileo/multibody/impulses/impulse-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct Impulse6dTpl;

    template <typename PhaseSpec>
    struct traits<Impulse6dTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = PS;

        using Meta_t = Impulse6dTpl<PS>;
        using Model_t = ImpulseModel6dTpl<PS>;
        using Data_t = ImpulseData6dTpl<PS>;

        using DimNC_t = DimensionTpl<6>;

        // Traits required by ForceDataBase
        using MatrixNcNv_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNV_t::Value, PS::Options>;
        using MatrixNcNdx_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, PS::DimNU_t::Value, PS::Options>;

        // Traits required by ImpulseDataBase
        using VectorNc_t = Eigen::GMatrix<typename PS::VarScalar, DimNC_t::Value, 1, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ImpulseData6dTpl<PhaseSpec>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = Impulse6dTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ImpulseModel6dTpl<PhaseSpec>>
    {
        using SpecOfBaseClass = PhaseSpec;
        using Meta_t = Impulse6dTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ImpulseData6dTpl
        : public ImpulseDataBase<ImpulseData6dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = Impulse6dTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ImpulseDataBase<ImpulseData6dTpl<PS>, PS>;

        GALILEO_IMPULSE_DATA_TYPEDEF(Meta_t);

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

        // Members required by ImpulseDataBase
        ActionMatrix_t fXj;
        MatrixNcNv_t dv0_dq;
        MatrixNv_t dtau_dq;

        // Accessor implementations required by ImpulseDataBase
        DEFAULT_ACCESSOR(ActionMatrix_t, fXj);
        DEFAULT_ACCESSOR(MatrixNcNv_t, dv0_dq);
        DEFAULT_ACCESSOR(MatrixNv_t, dtau_dq);

        // Members used for ImpulseModel6dTpl
        // Notice that we do not need to expose accessors for these because they are specific to the 6D impulse model
        SE3_t lwaMl;
        Motion_t v0;
        Force_t f_local;
        MatrixNcNv_t dv0_local_dq;
        Matrix6Nv_t fJf;
        Matrix6Nv_t v_partial_dq;
        Matrix6Nv_t v_partial_dv;
        Matrix3_t vv_skew;
        Matrix3_t vw_skew;
        Matrix3_t vv_world_skew;
        Matrix3_t vw_world_skew;
        Matrix3_t fv_skew;
        Matrix3_t fw_skew;
        MatrixNcNv_t fJf_df;

        ImpulseData6dTpl(const Model_t &model, RobotData_t *const robot_data)
            : robot(robot_data),
              frame(model.get_id()),
              type(model.get_type()),
              jMf(model.get_robot().frames[frame].placement),
              Jc(model.get_nc(), model.get_ps().get_nv()),
              f(Force_t::Zero()),
              fext(Force_t::Zero()),
              df_dx(model.get_nc(), model.get_ps().get_ndx()),
              df_du(model.get_nc(), model.get_ps().get_nu()),
              fXj(jMf.inverse().toActionMatrix()),
              dv0_dq(model.get_nc(), model.get_ps().get_nv()),
              dtau_dq(model.get_ps().get_nv(), model.get_ps().get_nv()),
              lwaMl(SE3_t::Identity()),
              v0(Motion_t::Zero()),
              f_local(Force_t::Zero()),
              dv0_local_dq(model.get_nc(), model.get_ps().get_nv()),
              fJf(6, model.get_ps().get_nv()),
              v_partial_dq(6, model.get_ps().get_nv()),
              v_partial_dv(6, model.get_ps().get_nv()),
              fJf_df(model.get_nc(), model.get_ps().get_nv())
        {
            Jc.setZero();
            df_dx.setZero();
            df_du.setZero();
            dv0_dq.setZero();
            dtau_dq.setZero();
            dv0_local_dq.setZero();
            fJf.setZero();
            v_partial_dq.setZero();
            v_partial_dv.setZero();
            vv_skew.setZero();
            vw_skew.setZero();
            vv_world_skew.setZero();
            vw_world_skew.setZero();
            fv_skew.setZero();
            fw_skew.setZero();
            fJf_df.setZero();
        }
    };

    template <typename PhaseSpec>
    struct ImpulseModel6dTpl
        : public ImpulseModelBase<ImpulseModel6dTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = Impulse6dTpl<PhaseSpec>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ImpulseModelBase<ImpulseModel6dTpl<PS>, PS>;

        using DimNC_t = typename traits<Meta_t>::DimNC_t;

        ImpulseModel6dTpl(const PS &ps,
                          const FrameIndex_t id,
                          const ReferenceFrame_t &type)
            : Base(ps, id, type, DimNC_t()),
              robot_(ps.get_state().get_robot())
        {
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            pinocchio::updateFramePlacement(get_robot(), *data.robot, get_id());
            pinocchio::getFrameJacobian(get_robot(), *data.robot, get_id(), pinocchio::LOCAL, data.fJf);

            switch (get_type())
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.Jc = data.fJf;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                data.lwaMl.rotation(data.robot->oMf[get_id()].rotation());
                data.Jc.noalias() = data.lwaMl.toActionMatrix() * data.fJf;
                break;
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const pinocchio::JointIndex joint = get_robot().frames[data.frame].parentJoint;
            pinocchio::getJointVelocityDerivatives(get_robot(), *data.robot,
                                                   joint, pinocchio::LOCAL,
                                                   data.v_partial_dq, data.v_partial_dv);
            data.dv0_local_dq.noalias() = data.fXj * data.v_partial_dq;

            switch (get_type())
            {
            case pinocchio::ReferenceFrame::LOCAL:
                data.dv0_dq = data.dv0_local_dq;
                break;
            case pinocchio::ReferenceFrame::WORLD:
            case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                const auto oRf = data.robot->oMf[get_id()].rotation();
                data.v0 = pinocchio::getFrameVelocity(get_robot(), *data.robot, get_id(),
                                                      pinocchio::LOCAL_WORLD_ALIGNED);
                pinocchio::skew(data.v0.linear(), data.vv_skew);
                pinocchio::skew(data.v0.angular(), data.vw_skew);
                data.vv_world_skew.noalias() = data.vv_skew * oRf;
                data.vw_world_skew.noalias() = data.vw_skew * oRf;
                data.dv0_dq.noalias() = data.lwaMl.toActionMatrix() * data.dv0_local_dq;
                topRows<3>(data.dv0_dq).noalias() -= data.vv_world_skew * bottomRows<3>(data.fJf);
                bottomRows<3>(data.dv0_dq).noalias() -= data.vw_world_skew * bottomRows<3>(data.fJf);
                break;
            }
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            data.f = Force_t(force);
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
                topRows<3>(data.fJf_df).noalias() = data.fv_skew * bottomRows<3>(data.fJf);
                bottomRows<3>(data.fJf_df).noalias() = data.fw_skew * bottomRows<3>(data.fJf);
                data.dtau_dq.noalias() = -data.fJf.transpose() * data.fJf_df;
                break;
            }
        }

        Data_t createData(RobotData_t *const robot) const
        {
            return Data_t(*this, robot);
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

    }; // struct ImpulseModel6dTpl

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_6d_hpp__
