#ifndef __galileo_predictive_segments_segment_erk_euler_hpp__
#define __galileo_predictive_segments_segment_erk_euler_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct SegmentERKEulerTpl;

    template <typename PhaseSpec>
    struct traits<SegmentERKEulerTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = SegmentERKModelEulerTpl<PS>;
        using Data_t = SegmentERKDataEulerTpl<PS>;

        using DimNStages_t = DimensionTpl<1>;
        static constexpr int NStages = DimNStages_t::Value;
    };

    template <typename PhaseSpec>
    struct traits<SegmentERKDataEulerTpl<PhaseSpec>>
    {
        using Meta_t = SegmentERKEulerTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<SegmentERKModelEulerTpl<PhaseSpec>>
    {
        using Meta_t = SegmentERKEulerTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct SegmentERKDataEulerTpl : public SegmentERKDataBase<SegmentERKDataEulerTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = SegmentERKDataBase<SegmentERKDataEulerTpl<PS>, PS>;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        DEFAULT_ACCESSOR(XNext_t, XNext);
        DEFAULT_ACCESSOR(XNextx_t, XNextx);
        DEFAULT_ACCESSOR(XNextw_t, XNextw);
        DEFAULT_ACCESSOR(L_t, L);
        DEFAULT_ACCESSOR(Lx_t, Lx);
        DEFAULT_ACCESSOR(Lw_t, Lw);
        DEFAULT_ACCESSOR(Lxx_t, Lxx);
        DEFAULT_ACCESSOR(Lxw_t, Lxw);
        DEFAULT_ACCESSOR(Lww_t, Lww);
        DEFAULT_ACCESSOR(H_t, H);
        DEFAULT_ACCESSOR(Hx_t, Hx);
        DEFAULT_ACCESSOR(Hw_t, Hw);
        DEFAULT_ACCESSOR(G_t, G);
        DEFAULT_ACCESSOR(Gx_t, Gx);
        DEFAULT_ACCESSOR(Gw_t, Gw);

        SegmentERKDataEulerTpl(const Model_t &model)
            : node(model.get_node().createData()),
              control(model.get_control().createData()),
              dx(model.get_ps().get_ndx()),
              da_dw(model.get_ps().get_nv(), model.get_ps().get_nw()),
              Luw(model.get_ps().get_nu(), model.get_ps().get_nw()),
              XNext(model.get_ps().get_nx()),
              XNextx(model.get_ps().get_ndx(), model.get_ps().get_ndx()),
              XNextw(model.get_ps().get_ndx(), model.get_ps().get_nw()),
              L(L_t(0.)),
              Lx(model.get_ps().get_ndx()),
              Lw(model.get_ps().get_nw()),
              Lxx(model.get_ps().get_ndx(), model.get_ps().get_ndx()),
              Lxw(model.get_ps().get_ndx(), model.get_ps().get_nw()),
              Lww(model.get_ps().get_nw(), model.get_ps().get_nw()),
              H(model.get_node().get_constraints().get_n_active()),
              Hx(model.get_node().get_constraints().get_n_active(), model.get_ps().get_ndx()),
              Hw(model.get_node().get_constraints().get_n_active(), model.get_ps().get_nw()),
              G(model.get_node().get_constraints().get_n_active()),
              Gx(model.get_node().get_constraints().get_n_active(), model.get_ps().get_ndx()),
              Gw(model.get_node().get_constraints().get_n_active(), model.get_ps().get_nw())
        {
            dx.setZero();
            da_dw.setZero();
            Luw.setZero();

            XNext.setZero();
            XNextx.setZero();
            XNextw.setZero();
            Lx.setZero();
            Lw.setZero();
            Lxx.setZero();
            Lxw.setZero();
            Lww.setZero();
            H.setZero();
            Hx.setZero();
            Hw.setZero();
            G.setZero();
            Gx.setZero();
            Gw.setZero();
        }

        NodeData_t node;
        ControlParamData_t control;
        VectorNdx_t dx;
        MatrixNvNw_t da_dw;
        MatrixNuNw_t Luw;

        XNext_t XNext;
        XNextx_t XNextx;
        XNextw_t XNextw;

        L_t L;
        Lx_t Lx;
        Lw_t Lw;
        Lxx_t Lxx;
        Lxw_t Lxw;
        Lww_t Lww;

        H_t H;
        Hx_t Hx;
        Hw_t Hw;

        G_t G;
        Gx_t Gx;
        Gw_t Gw;

    }; // struct SegmentERKDataEulerTpl

    template <typename PhaseSpec>
    class SegmentERKModelEulerTpl : public SegmentERKModelBase<SegmentERKModelEulerTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = SegmentERKModelBase<SegmentERKModelEulerTpl<PS>, PS>;

        SegmentERKModelEulerTpl(const PS &ps,
                                const NodeModel_t &node,
                                const ControlParamModel_t &control,
                                const NumScalar period)
            : Base(ps), node_(node), control_(control), period_(period), period_squared_(period * period)
        {
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            const auto nv_dim = get_ps().get_nv_dim();
            const auto v = tail(x, nv_dim);
            control_.calc(data.control, NumScalar(0.), w);
            node_.calc(data.node, x, data.control.u);
            const VectorNv_t &a = data.node.XAcc_accessor();
            head(data.dx, nv_dim).noalias() = v * period_ + a * period_squared_;
            tail(data.dx, nv_dim).noalias() = a * period_;

            get_state().integrate(x, data.dx, data.XNext);
            data.L = period_ * data.node.L_accessor();
            data.H = data.node.H_accessor();
            data.G = data.node.G_accessor();
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            node_.calc(data.node, x);
            data.dx.setZero();
            data.XNext = x;
            data.L = data.node.L_accessor();
            data.H = data.node.H_accessor();
            data.G = data.node.G_accessor();
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            control_.calc(data.control, NumScalar(0.), w);
            node_.calcDiff(data.node, x, data.control.u);
            const MatrixNvNdx_t &da_dx = data.node.XAccx_accessor();
            const MatrixNvNw_t &da_du = data.node.XAccu_accessor();
            control_.multiplyByJacobian(data.control, da_du, data.da_dw);

            const auto nv_dim = get_ps().get_nv_dim();
            topRows(data.XNextx, nv_dim).noalias() = da_dx * period_squared_;
            bottomRows(data.XNextx, nv_dim).noalias() = da_dx * period_;
            topRightCorner(data.XNextx, nv_dim, nv_dim).diagonal().array() += VarScalar(period_);
            topRows(data.XNextw, nv_dim).noalias() = period_squared_ * data.da_dw;
            bottomRows(data.XNextw, nv_dim).noalias() = period_ * data.da_dw;

            const auto state = get_state();
            state.template JintegrateTransport<SECOND>(x, data.dx, data.XNextx);
            state.template Jintegrate<FIRST, ADDTO>(x, data.dx, data.XNextx, data.XNextx);
            state.template JintegrateTransport<SECOND>(x, data.dx, data.XNextw);

            data.Lx.noalias() = period_ * data.node.Lx_accessor();
            control_.multiplyJacobianTransposeBy(data.control, data.node.Lu_accessor(), data.Lw);
            data.Lw *= period_;
            data.Lxx.noalias() = period_ * data.node.Lxx_accessor();
            control_.multiplyByJacobian(data.control, data.node.Lxu_accessor(), data.Lxw);
            data.Lxw *= period_;
            control_.multiplyByJacobian(data.control, data.node.Luu_accessor(), data.Luw);
            control_.multiplyJacobianTransposeBy(data.control, data.Luw, data.Lww);
            data.Lww *= period_;

            data.Gx = data.node.Gx_accessor();
            data.Hx = data.node.Hx_accessor();
            data.Gw.conservativeResize(node_.get_constraints().get_n_active(), get_ps().get_nw());
            data.Hw.conservativeResize(node_.get_constraints().get_n_active(), get_ps().get_nw());
            control_.multiplyByJacobian(data.control, data.node.Gu_accessor(), data.Gw);
            control_.multiplyByJacobian(data.control, data.node.Hu_accessor(), data.Hw);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            node_.calcDiff(data.node, x);
            get_state().Jintegrate(x, data.dx, data.XNextx, data.XNextx);
            data.Lx = data.node.Lx_accessor();
            data.Lxx = data.node.Lxx_accessor();
            data.Gx = data.node.Gx_accessor();
            data.Hx = data.node.Hx_accessor();
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlParamVectorType> &w,
                         const int maxiter,
                         const NumScalar tol) const
        {
            data.control.u.setZero();
            node_.quasiStatic(data.node, x, data.control.u, maxiter, tol);
            control_.params(data.control, NumScalar(0.), data.control.u);
            w = data.control.w;
        }

        Data_t createData() const { return Data_t(*this); }

        const ControlParamModel_t &get_control() const { return control_; }
        const NodeModel_t &get_node() const { return node_; }

        using Base::get_ps;
        using Base::get_state;

    protected:
        NodeModel_t node_;
        ControlParamModel_t control_;

        NumScalar period_;
        NumScalar period_squared_;

    }; // class SegmentERKModelEulerTpl

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_euler_hpp__
