#ifndef __galileo_core_controls_control_param_polynomial_hpp__
#define __galileo_core_controls_control_param_polynomial_hpp__

#include "galileo/core/controls/control-param-base.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"

namespace galileo
{

    template <typename PhaseSpec, int NOrder_>
    struct ControlParamPolynomialTpl;

    template <typename PhaseSpec, int NOrder_>
    struct traits<ControlParamPolynomialTpl<PhaseSpec, NOrder_>>
    {
        using PS = PhaseSpec;

        using Meta_t = ControlParamPolynomialTpl<PS, NOrder_>;
        using Model_t = ControlParamModelPolynomialTpl<PS, NOrder_>;
        using Data_t = ControlParamDataPolynomialTpl<PS, NOrder_>;

        using DimNOrder_t = DimensionTpl<NOrder_>;
        static constexpr int NOrder = DimNOrder_t::Value;
    };

    template <typename PhaseSpec, int NOrder_>
    struct traits<ControlParamModelPolynomialTpl<PhaseSpec, NOrder_>>
    {
        using Meta_t = ControlParamPolynomialTpl<PhaseSpec, NOrder_>;
    };

    template <typename PhaseSpec, int NOrder_>
    struct traits<ControlParamDataPolynomialTpl<PhaseSpec, NOrder_>>
    {
        using Meta_t = ControlParamPolynomialTpl<PhaseSpec, NOrder_>;
    };

    template <typename PhaseSpec, int NOrder_>
    class ControlParamDataPolynomialTpl
        : public ControlParamDataBase<ControlParamDataPolynomialTpl<PhaseSpec, NOrder_>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ControlParamPolynomialTpl<PS, NOrder_>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ControlParamDataBase<ControlParamDataPolynomialTpl<PS, NOrder_>, PS>;

        using VectorNu_t = ArenaMatrixTpl<typename PS::VectorNu_t>;
        using VectorNw_t = ArenaMatrixTpl<typename PS::VectorNw_t>;
        using MatrixNuNw_t = ArenaMatrixTpl<typename PS::MatrixNuNw_t>;

        DEFAULT_ACCESSOR(VectorNu_t, u);
        DEFAULT_ACCESSOR(VectorNw_t, w);
        DEFAULT_ACCESSOR(MatrixNuNw_t, du_dw);

        ControlParamDataPolynomialTpl(const Model_t &model, MemoryArena &arena)
            : u(arena, model.get_ps().get_nu(), 1),
              w(arena, model.get_ps().get_nw(), 1),
              du_dw(arena, model.get_ps().get_nu(), model.get_ps().get_nw())
        {
            u.setZero();
            w.setZero();
            du_dw.setZero();
        }

        VectorNu_t u;
        VectorNw_t w;
        MatrixNuNw_t du_dw;

    }; // class ControlParamDataPolynomialTpl

    template <typename PhaseSpec, int NOrder_>
    class ControlParamModelPolynomialTpl
        : public ControlParamModelBase<ControlParamModelPolynomialTpl<PhaseSpec, NOrder_>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ControlParamPolynomialTpl<PS, NOrder_>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ControlParamModelBase<ControlParamModelPolynomialTpl<PS, NOrder_>, PS>;

        using NumScalar = typename PS::NumScalar;
        using BarycentricInterpolator_t = BarycentricInterpolatorTpl<NumScalar, NOrder_, PS::Options>;

        ControlParamModelPolynomialTpl(const PS &ps, const BarycentricInterpolator_t &interpolator)
            : Base(ps), interpolator_(interpolator)
        {
        }

        template <typename ControlParamVectorType>
        void calc(Data_t &data, const NumScalar &t, const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            interpolator_.calc(t, w.reshaped(get_ps().get_nu(), get_ps().get_norder()), data.u);
        }

        template <typename ControlParamVectorType>
        void calcDiff(Data_t &data, const NumScalar &t, const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            interpolator_.calcDiff(t, w.reshaped(get_ps().get_nu(), get_ps().get_norder()), data.du_dw);
        }

        template <typename ControlVectorType>
        void params(Data_t &data, const NumScalar &t, const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            for (int i = 0; i < get_ps().get_norder(); ++i)
                segment(data.w, i * get_ps().get_nu(), get_ps().get_nu_dim()) = u;
        }

        template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
        void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                           const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
        {
            for (int i = 0; i < get_ps().get_norder(); ++i)
            {
                segment(w_lb, i * get_ps().get_nu(), get_ps().get_nu_dim()) = u_lb;
                segment(w_ub, i * get_ps().get_nu(), get_ps().get_nu_dim()) = u_ub;
            }
        }

        template <AssignmentOp op = SETTO, typename InputMatrixType, typename OutputMatrixType>
        void multiplyByJacobian(Data_t &data,
                                const Eigen::MatrixBase<InputMatrixType> &A,
                                Eigen::MatrixBase<OutputMatrixType> &out) const
        {
            DimensionTpl<InputMatrixType::RowsAtCompileTime> A_rows_dim(A.rows());
            for (int i = 0; i < get_ps().get_norder(); ++i)
            {
                auto out_block_i = block(out, 0, i * get_ps().get_nu(), A_rows_dim, get_ps().get_nu_dim());
                const auto du_dw_block_i =
                    block(data.du_dw, 0, i * get_ps().get_nu(), get_ps().get_nu_dim(), get_ps().get_nu_dim());
                if constexpr (IsSetTo<op>)
                    out_block_i.noalias() = A * du_dw_block_i;
                else if constexpr (IsAddTo<op>)
                    out_block_i += A * du_dw_block_i;
                else if constexpr (IsRmFrom<op>)
                    out_block_i -= A * du_dw_block_i;
            }
        }

        template <AssignmentOp op = SETTO, typename InputMatrixType, typename OutputMatrixType>
        void multiplyJacobianTransposeBy(Data_t &data,
                                         const Eigen::MatrixBase<InputMatrixType> &A,
                                         Eigen::MatrixBase<OutputMatrixType> &out) const
        {
            DimensionTpl<InputMatrixType::ColsAtCompileTime> A_cols_dim(A.cols());
            for (int i = 0; i < get_ps().get_norder(); ++i)
            {
                auto out_block_i = block(out, i * get_ps().get_nu(), 0, get_ps().get_nu_dim(), A_cols_dim);
                const auto du_dw_block_i =
                    block(data.du_dw, 0, i * get_ps().get_nu(), get_ps().get_nu_dim(), get_ps().get_nu_dim());
                if constexpr (IsSetTo<op>)
                    out_block_i.noalias() = du_dw_block_i.transpose() * A;
                else if constexpr (IsAddTo<op>)
                    out_block_i += du_dw_block_i.transpose() * A;
                else if constexpr (IsRmFrom<op>)
                    out_block_i -= du_dw_block_i.transpose() * A;
            }
        }

        Data_t createData(MemoryArena &arena) const { return Data_t(*this, arena); }

        using Base::get_ps;

    protected:
        BarycentricInterpolator_t interpolator_;

    }; // class ControlParamModelPolynomialTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_polynomial_hpp__
