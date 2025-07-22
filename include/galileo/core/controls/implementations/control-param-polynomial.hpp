#ifndef __galileo_core_controls_control_param_polynomial_hpp__
#define __galileo_core_controls_control_param_polynomial_hpp__

#include "galileo/core/controls/control-param-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

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
        using Data_t = ControlParamDataTpl<PS>;

        using DimNOrder_t = DimensionTpl<NOrder_>;
    };

    template <typename PhaseSpec, int NOrder_>
    struct traits<ControlParamModelPolynomialTpl<PhaseSpec, NOrder_>>
    {
        using PS = PhaseSpec;

        using Meta_t = ControlParamPolynomialTpl<PS, NOrder_>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec, int NOrder_>
    class ControlParamModelPolynomialTpl
        : public ControlParamModelBase<ControlParamModelPolynomialTpl<PhaseSpec, NOrder_>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ControlParamPolynomialTpl<PS, NOrder_>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ControlParamModelBase<ControlParamModelPolynomialTpl<PS, NOrder_>, PS>;

        GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PS);
        GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PS);

        using BarycentricInterpolator_t = BarycentricInterpolatorTpl<NumScalar, NOrder_, Options>;

        ControlParamModelPolynomialTpl(const PS &ps, const BarycentricInterpolator_t &interpolator)
            : Base(ps),
              interpolator_(interpolator)
        {
        }

        template <typename ControlParamVectorType>
        void calc(Data_t &data, const NumScalar &t,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            interpolator_.calc(t, w.reshaped(get_nu(), get_norder()), data.u);
        }

        template <typename ControlParamVectorType>
        void calcDiff(Data_t &data, const NumScalar &t,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            interpolator_.calcDiff(t, w.reshaped(get_nu(), get_norder()), data.du_dw);
        }

        template <typename ControlVectorType>
        void params(Data_t &data, const NumScalar &t,
                    const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            for (int i = 0; i < get_norder(); ++i)
            {
                segment(data.w, i * get_nu(), get_nu_dim()) = u;
            }
        }

        template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
        void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                           const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
        {
            for (int i = 0; i < get_norder(); ++i)
            {
                segment(w_lb, i * get_nu(), get_nu_dim()) = u_lb;
                segment(w_ub, i * get_nu(), get_nu_dim()) = u_ub;
            }
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyByJacobian(
            Data_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            switch (op)
            {
            case setto:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, 0, i * get_nu(), get_nw_dim(), get_nu_dim()) = data.du_dw(0, i * get_nu()) * A;
                }
                break;
            case addto:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, 0, i * get_nu(), get_nw_dim(), get_nu_dim()) += data.du_dw(0, i * get_nu()) * A;
                }
                break;
            case rmfrom:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, 0, i * get_nu(), get_nw_dim(), get_nu_dim()) -= data.du_dw(0, i * get_nu()) * A;
                }
                break;
            default:
                break;
            }
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyJacobianTransposeBy(
            Data_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            switch (op)
            {
            case setto:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, i * get_nu(), 0, get_nu_dim(), get_nw_dim()) = data.du_dw(0, i * get_nu()) * A;
                }
                break;
            case addto:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, i * get_nu(), 0, get_nu_dim(), get_nw_dim()) += data.du_dw(0, i * get_nu()) * A;
                }
                break;
            case rmfrom:
                for (int i = 0; i < get_norder(); ++i)
                {
                    block(out, i * get_nu(), 0, get_nu_dim(), get_nw_dim()) -= data.du_dw(0, i * get_nu()) * A;
                }
                break;
            default:
                break;
            }
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        using Base::get_ps;

        using Base::get_nu;
        using Base::get_nu_dim;

        using Base::get_norder;
        using Base::get_norder_dim;

        using Base::get_nw;
        using Base::get_nw_dim;

    protected:
        BarycentricInterpolatorTpl<NumScalar, PS::NOrder, PS::Options> interpolator_;

    }; // class ControlParamModelPolynomialTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_polynomial_hpp__
