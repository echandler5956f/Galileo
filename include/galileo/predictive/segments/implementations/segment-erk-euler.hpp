#ifndef __galileo_predictive_segments_segment_erk_euler_hpp__
#define __galileo_predictive_segments_segment_erk_euler_hpp__

#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{

    template <
        typename PhaseSpec>
    struct SegmentERKEulerTpl;

    template <
        typename PhaseSpec>
    struct traits<SegmentERKEulerTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = SegmentERKDataEulerTpl<PS>;
        using Data_t = SegmentERKModelEulerTpl<PS>;

        static constexpr int NStages = 1;
    };

    template <
        typename PhaseSpec>
    struct traits<SegmentERKDataEulerTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <
        typename PhaseSpec>
    struct traits<SegmentERKModelEulerTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct SegmentERKDataEulerTpl : public SegmentERKDataBase<SegmentERKDataEulerTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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

        NodeData_t node;
        ControlParamData_t control;
        VectorNdx_t dx;
        MatrixNvNw_t da_dw;

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
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = SegmentERKEulerTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        template <typename StateVectorType, typename ControlParamVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NV> v =
                x.tail(NV);

            control_.calc(data.control, NumScalar(0.), w);
            node_.calc(data.node, x, data.control.u);
            const VectorNdx_t &a = data.node.Xacc;
            data.dx.head(NV).noalias() = v * period_ + a * period_squared_;
            data.dx.tail(NV).noalias() = a * period_;
            state_.integrate(x, data.dx, data.XNext);
            data.L = period_ * data.node.L;
            data.G = data.node.G;
            data.H = data.node.H;
            if (with_cost_residual_)
            {
                data.r = data.node.r;
            }
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            node_.calc(data.node, x);
            data.dx.setZero();
            data.XNext = x;
            data.L = data.node.L;
            data.G = data.node.G;
            data.H = data.node.H;
            if (with_cost_residual_)
            {
                data.r = data.node.r;
            }
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            control_.calc(data.control, NumScalar(0.), w);
            node_.calcDiff(data.node, x, data.control.u);
            const MatrixNvNdx_t &da_dx = data.node.Fx;
            const MatrixNvNw_t &da_dw = data.node.Fu;
            control_.multiplyByJacobian(data.control, da_dw, data.da_dw);
            data.XNextx.topRows(NV).noalias() = da_dx * period_squared_;
            data.XNextx.bottomRows(NV).noalias() = da_dx * period_;
            data.XNextx.topRightCorner(NV, NV).diagonal().array() += NumScalar(period_);
            data.XNextw.topRows(NV).noalias() = period_squared_ * data.da_dw;
            data.XNextw.bottomRows(NV).noalias() = period_ * data.da_dw;
            state_.JintegrateTransport(x, data.dx, data.XNextx, second);
            state_.Jintegrate(x, data.dx, data.XNextx, data.XNextx, first, addto);
            state_.JintegrateTransport(x, data.dx, data.XNextw, second);

            data.Lx.noalias() = period_ * data.node.Lx;
            control_.multiplyJacobianTransposeBy(data.control, data.node.Lu, data.Lw);
            data.Lw *= period_;
            data.Lxx.noalias() = period_ * data.node.Lxx;
            control_.multiplyByJacobian(data.control, data.node.Lxu, data.Lxw);
            data.Lxw *= period_;
            control_.multiplyByJacobian(data.control, data.node.Luu, data.Luw);
            control_.multiplyJacobianTransposeBy(data.control, data.Luw, data.Lww);
            data.Lww *= period_;
            data.Gx = data.node.Gx;
            data.Hx = data.node.Hx;
            data.Gw.resize(node_.ng(), NW);
            data.Hw.resize(node_.nh(), NW);
            control_.multiplyByJacobian(data.control, data.node.Gu, data.Gw);
            control_.multiplyByJacobian(data.control, data.node.Hu, data.Hw);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            node_.calcDiff(data.node, x);
            state_.Jintegrate(x, data.dx, data.XNextx, data.XNextx);
            data.Lx = data.node.Lx;
            data.Lxx = data.node.Lxx;
            data.Gx = data.node.Gx;
            data.Hx = data.node.Hx;
        }

        template <typename StateVectorType, typename ControlParamVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlParamVectorType> &w,
                         const int maxiter, const NumScalar tol) const
        {
            data.control.u.setZero();
            node_.quasiStatic(data.node, x, data.control.u, maxiter, tol);
            control_.params(data.control, NumScalar(0.), data.control.u);
            w = data.control.w;
        }

    protected:
        State_t state_;
        ControlParamModel_t control_;
        NodeModel_t node_;

        NumScalar period_;
        NumScalar period_squared_;

        bool with_cost_residual_;

    }; // class SegmentERKModelEulerTpl

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_euler_hpp__
