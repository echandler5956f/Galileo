#ifndef __galileo_multibody_predictive_jumps_jump_impulse_fwddyn_hpp__
#define __galileo_multibody_predictive_jumps_jump_impulse_fwddyn_hpp__

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/impulse-dynamics.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include "galileo/core/data/data-collector-default.hpp"
#include "galileo/domains/multibody/spatial/impulses/impulse-manager.hpp"
#include "galileo/predictive/jumps/jump-base.hpp"
#include "galileo/domains/multibody/predictive/jumps/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct JumpImpulseFwdDynTpl;

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<JumpImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = JumpImpulseFwdDynTpl<PS, ImpulseCollectionTpl>;
        using Model_t = JumpModelImpulseFwdDynTpl<PS, ImpulseCollectionTpl>;
        using Data_t = JumpDataImpulseFwdDynTpl<PS, ImpulseCollectionTpl>;

        using DimNC_t = DimensionTpl<Eigen::Dynamic>;
        static constexpr int NC = DimNC_t::Value;

        using ImpulseManagerMeta_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using ImpulseModelManager_t = typename traits<ImpulseManagerMeta_t>::ModelManager_t;
        using ImpulseDataManager_t = typename traits<ImpulseManagerMeta_t>::DataManager_t;

        using DataCollector_t = DataCollectorImpulseTpl<PS, ImpulseCollectionTpl>;

        using Kinv_t =
            ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, AddDim_v<PS::NV, NC>, AddDim_v<PS::NV, NC>>>;
        using MatrixNcNdx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NC, PS::NDX>>;
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<JumpDataImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using Meta_t = JumpImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>;
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct traits<JumpModelImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using Meta_t = JumpImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>;
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct JumpDataImpulseFwdDynTpl
        : public JumpDataBase<JumpDataImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = JumpImpulseFwdDynTpl<PS, ImpulseCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = JumpDataBase<JumpDataImpulseFwdDynTpl<PS, ImpulseCollectionTpl>, PS>;

        using ImpulseManagerMeta_t = typename traits<Meta_t>::ImpulseManagerMeta_t;
        using ImpulseModelManager_t = typename traits<Meta_t>::ImpulseModelManager_t;
        using ImpulseDataManager_t = typename traits<Meta_t>::ImpulseDataManager_t;

        using XNext_t = ArenaMatrixTpl<typename PS::XNext_t>;
        using XNextx_t = ArenaMatrixTpl<typename PS::XNextx_t>;
        using L_t = typename PS::L_t;
        using Lx_t = ArenaMatrixTpl<typename PS::Lx_t>;
        using Lxx_t = ArenaMatrixTpl<typename PS::Lxx_t>;
        using H_t = ArenaMatrixTpl<typename PS::H_t>;
        using Hx_t = ArenaMatrixTpl<typename PS::Hx_t>;
        using G_t = ArenaMatrixTpl<typename PS::G_t>;
        using Gx_t = ArenaMatrixTpl<typename PS::Gx_t>;

        using DataCollector_t = typename traits<Meta_t>::DataCollector_t;
        using RobotData_t = typename DataCollector_t::RobotData_t;
        using CostDataManager_t = typename PS::CostDataManager_t;
        using ConstraintDataManager_t = typename PS::ConstraintDataManager_t;

        using VectorNv_t = ArenaMatrixTpl<typename PS::VectorNv_t>;
        using MatrixNv_t = ArenaMatrixTpl<typename PS::MatrixNv_t>;

        using Kinv_t = typename traits<Meta_t>::Kinv_t;
        using MatrixNcNdx_t = typename traits<Meta_t>::MatrixNcNdx_t;

        DEFAULT_ACCESSOR(XNext_t, XNext);
        DEFAULT_ACCESSOR(XNextx_t, XNextx);

        L_t &L_accessor() { return costs.L; }
        const L_t &L_accessor() const { return costs.L; }
        Lx_t &Lx_accessor() { return costs.Lx; }
        const Lx_t &Lx_accessor() const { return costs.Lx; }
        Lxx_t &Lxx_accessor() { return costs.Lxx; }
        const Lxx_t &Lxx_accessor() const { return costs.Lxx; }
        H_t &H_accessor() { return constraints.H; }
        const H_t &H_accessor() const { return constraints.H; }
        Hx_t &Hx_accessor() { return constraints.Hx; }
        const Hx_t &Hx_accessor() const { return constraints.Hx; }
        G_t &G_accessor() { return constraints.H; }
        const G_t &G_accessor() const { return constraints.H; }
        Gx_t &Gx_accessor() { return constraints.Hx; }
        const Gx_t &Gx_accessor() const { return constraints.Hx; }

        JumpDataImpulseFwdDynTpl(const Model_t &model, MemoryArena &arena)
            : XNext(arena, model.get_ps().get_nx()),
              XNextx(arena, model.get_ps().get_ndx(), model.get_ps().get_ndx()),
              data_collector(
                  std::make_shared<DataCollector_t>(std::make_shared<RobotData_t>(model.get_state().get_robot()))),
              robot(data_collector->robot),
              impulses(model.get_impulses().createData(arena, robot.get())),
              costs(model.get_costs().createData(arena, data_collector.get())),
              constraints(model.get_constraints().createData(arena, data_collector.get())),
              vnone(arena, model.get_ps().get_nv()),
              Kinv(arena,
                   model.get_ps().get_nv() + model.get_impulses().get_n_total(),
                   model.get_ps().get_nv() + model.get_impulses().get_n_total()),
              df_dx(arena, model.get_impulses().get_n_total(), model.get_ps().get_ndx()),
              dgrav_dq(arena, model.get_ps().get_nv(), model.get_ps().get_nv())
        {
            XNext.setZero();
            XNextx.setZero();
            vnone.setZero();
            Kinv.setZero();
            df_dx.setZero();
            dgrav_dq.setZero();
        }

        XNext_t XNext;
        XNextx_t XNextx;

        std::shared_ptr<DataCollector_t> data_collector;
        std::shared_ptr<RobotData_t> robot;
        ImpulseDataManager_t impulses;
        CostDataManager_t costs;
        ConstraintDataManager_t constraints;

        VectorNv_t vnone;
        Kinv_t Kinv;
        MatrixNcNdx_t df_dx;
        MatrixNv_t dgrav_dq;

    }; // class JumpDataImpulseFwdDynTpl

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    class JumpModelImpulseFwdDynTpl
        : public JumpModelBase<JumpModelImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = JumpImpulseFwdDynTpl<PS, ImpulseCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = JumpModelBase<JumpModelImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTpl>, PhaseSpec>;

        using ImpulseManagerMeta_t = typename traits<Meta_t>::ImpulseManagerMeta_t;
        using ImpulseModelManager_t = typename traits<Meta_t>::ImpulseModelManager_t;
        using ImpulseDataManager_t = typename traits<Meta_t>::ImpulseDataManager_t;

        using CostModelManager_t = typename PS::CostModelManager_t;
        using ConstraintModelManager_t = typename PS::ConstraintModelManager_t;
        using State_t = typename PS::State_t;
        using NumScalar = typename PS::NumScalar;
        using VectorNv_t = typename PS::VectorNv_t;
        using Motion_t = typename PS::Motion_t;

        JumpModelImpulseFwdDynTpl(PS &ps,
                                  const State_t &state,
                                  const CostModelManager_t &costs,
                                  const ConstraintModelManager_t &constraints,
                                  const ImpulseModelManager_t &impulses,
                                  const NumScalar &r_coeff = 0.0,
                                  const NumScalar &JMinvJt_damping = 0.0,
                                  const bool enable_force = false)
            : Base(ps, state),
              costs_(costs),
              constraints_(constraints),
              impulses_(impulses),
              with_armature_(true),
              armature_(get_ps().get_nv()),
              r_coeff_(fabs(r_coeff)),
              JMinvJt_damping_(fabs(JMinvJt_damping)),
              enable_force_(enable_force),
              gravity_(get_state().get_robot().gravity)
        {
            armature_.setZero();
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const auto nc_dim = get_impulses().get_n_active_dim();

            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());
            pinocchio::computeAllTerms(get_state().get_robot(), *data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_state().get_robot(), *data.robot);
            if (!with_armature_)
            {
                data.robot->M.diagonal() += armature_;
            }
            get_impulses().calc(data.impulses, x);

            pinocchio::impulseDynamics(get_state().get_robot(),
                                       *data.robot,
                                       v,
                                       topRows(data.impulses.Jc, nc_dim),
                                       r_coeff_,
                                       JMinvJt_damping_);

            head(data.XNext, get_ps().get_nq_dim()) = q;
            tail(data.XNext, get_ps().get_nv_dim()) = data.robot->dq_after;
            get_impulses().updateVelocity(data.impulses, data.robot->dq_after);
            get_impulses().updateForce(data.impulses, data.robot->impulse_c);

            get_costs().calc(data.costs, x);
            data.L_accessor() = data.costs.L;
            if (get_constraints().get_n_active() > 0)
            {
                get_constraints().calc(data.constraints, x);
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const auto nv_dim = get_ps().get_nv_dim();
            const auto nc_dim = get_impulses().get_n_active_dim();

            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());

            // Computing the dynamics derivatives
            // We resize the Kinv matrix because Eigen cannot call block operations
            // recursively: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=408. Therefore,
            // it is not possible to pass data.Kinv.topLeftCorner(nv + nc, nv + nc)
            data.Kinv.resize(get_ps().get_nv() + nc_dim.value(), nv_dim.value() + nc_dim.value());
            pinocchio::computeRNEADerivatives(get_state().get_robot(),
                                              *data.robot,
                                              q,
                                              data.vnone,
                                              data.robot->dq_after - v,
                                              data.impulses.fext);
            pinocchio::computeGeneralizedGravityDerivatives(
                get_state().get_robot(), *data.robot, q, data.dgrav_dq);
            pinocchio::getKKTContactDynamicMatrixInverse(
                get_state().get_robot(), *data.robot, topRows(data.impulses.Jc, nc_dim), data.Kinv);

            pinocchio::computeForwardKinematicsDerivatives(
                get_state().get_robot(), *data.robot, q, data.robot->dq_after, data.vnone);
            get_impulses().calcDiff(data.impulses, x);
            get_impulses().updateRneaDiff(data.impulses, *data.robot);

            const auto a_partial_dtau = topLeftCorner(data.Kinv, nv_dim, nv_dim);
            const auto a_partial_da = topRightCorner(data.Kinv, nv_dim, nc_dim);
            const auto f_partial_dtau = bottomLeftCorner(data.Kinv, nc_dim, nv_dim);
            const auto f_partial_da = bottomRightCorner(data.Kinv, nc_dim, nc_dim);

            data.robot->dtau_dq -= data.dgrav_dq;
            data.robot->M.template triangularView<Eigen::StrictlyLower>() =
                data.robot->M.transpose().template triangularView<Eigen::StrictlyLower>();
            topLeftCorner(data.XNextx, nv_dim, nv_dim).setIdentity();
            topRightCorner(data.XNextx, nv_dim, nv_dim).setZero();
            bottomLeftCorner(data.XNextx, nv_dim, nv_dim).noalias() = -a_partial_dtau * data.robot->dtau_dq;
            bottomLeftCorner(data.XNextx, nv_dim, nv_dim).noalias() -=
                a_partial_da * topRows(data.impulses.dv0_dq, nc_dim);
            bottomRightCorner(data.XNextx, nv_dim, nv_dim).noalias() = a_partial_dtau * data.robot->M;

            // Computing the cost derivatives
            if (enable_force_)
            {
                topLeftCorner(data.df_dx, nc_dim, nv_dim).noalias() = f_partial_dtau * data.robot->dtau_dq;
                topLeftCorner(data.df_dx, nc_dim, nv_dim).noalias() +=
                    f_partial_da * topRows(data.impulses.dv0_dq, nc_dim);
                topRightCorner(data.df_dx, nc_dim, nv_dim).noalias() = f_partial_da * topRows(data.impulses.Jc, nc_dim);
                get_impulses().updateVelocityDiff(data.impulses, bottomRows(data.XNextx, nv_dim));
                get_impulses().updateForceDiff(data.impulses, topRows(data.df_dx, nc_dim));
            }

            get_costs().calcDiff(data.costs, x);
            if (get_constraints().get_n_active() > 0 || get_constraints().get_n_active() > 0)
            {
                get_constraints().calcDiff(data.constraints, x);
            }
        }

        Data_t createData(MemoryArena &arena) const { return Data_t(*this, arena); }

        const CostModelManager_t &get_costs() const { return costs_.get(); }
        const ConstraintModelManager_t &get_constraints() const { return constraints_.get(); }
        const ImpulseModelManager_t &get_impulses() const { return impulses_.get(); }

        using Base::get_ps;
        using Base::get_state;

    protected:
        std::reference_wrapper<const CostModelManager_t> costs_;
        std::reference_wrapper<const ConstraintModelManager_t> constraints_;
        std::reference_wrapper<const ImpulseModelManager_t> impulses_;

        bool with_armature_ = true;
        VectorNv_t armature_;
        NumScalar r_coeff_;
        NumScalar JMinvJt_damping_;
        bool enable_force_ = false;
        Motion_t gravity_;

    }; // class JumpModelImpulseFwdDynTpl

} // namespace galileo

#endif // __galileo_multibody_predictive_jumps_jump_impulse_fwddyn_hpp__
