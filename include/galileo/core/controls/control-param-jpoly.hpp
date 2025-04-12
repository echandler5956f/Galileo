#ifndef __galileo_core_controls_control_param_jpoly_hpp__
#define __galileo_core_controls_control_param_jpoly_hpp__

#include "galileo/core/controls/controls-param-base.hpp"
#include "galileo/math/polynomial.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ControlParamDataJacobiPolynomialTpl : public ControlParamDataBase<ControlParamDataJacobiPolynomialTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        typename PS::U_t u;
        typename PS::W_t w;
        typename PS::Uw_t du_dw;

    }; // struct ControlParamDataJacobiPolynomialTpl

    template <typename PhaseSpec>
    class ControlParamModelJacobiPolynomialTpl : public ControlParamModelBase<ControlParamModelJacobiPolynomialTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        template <typename ControlParamVectorType>
        void calc(typename PS::ControlParamData_t &data, const typename PS::NumScalar &t,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            jacobi_polynomial_.barycentricInterpolation(t, w.reshaped(PS::NU, PS::NOrder), data.u.derived());
        }

        template <typename ControlParamVectorType>
        void calcDiff(typename PS::ControlParamData_t &data,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            jacobi_polynomial_.barycentricInterpolationDiff(t, w.reshaped(PS::NU, PS::NOrder), data.du_dw.derived());
        }

        template <typename ControlVectorType>
        void params(typename PS::ControlParamData_t &data, const typename PS::NumScalar &t,
                    const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            for (std::size_t i = 0; i < PS::NOrder; ++i)
            {
                data.w.segment(i * PS::NU, PS::NU) = u;
            }
        }

        template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
        void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                           const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
        {
            for (std::size_t i = 0; i < PS::NOrder; ++i)
            {
                w_lb.segment(i * PS::NU, PS::NU) = u_lb;
                w_ub.segment(i * PS::NU, PS::NU) = u_ub;
            }
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyByJacobian(
            typename PS::ControlParamData_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            switch (op)
            {
            case setto:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(0, i * PS::NU, PS::NW, PS::NU) = typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            case addto:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(0, i * PS::NU, PS::NW, PS::NU) += typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            case rmfrom:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(0, i * PS::NU, PS::NW, PS::NU) -= typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            default:
                break;
            }
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyJacobianTransposeBy(
            typename PS::ControlParamData_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            switch (op)
            {
            case setto:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(i * PS::NU, 0, PS::NU, PS::NW) = typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            case addto:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(i * PS::NU, 0, PS::NU, PS::NW) += typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            case rmfrom:
                for (std::size_t i = 0; i < PS::NOrder; ++i)
                {
                    out.block(i * PS::NU, 0, PS::NU, PS::NW) -= typename PS::NumScalar(data.du_dw(0, i * PS::NU)) * A;
                }
                break;
            default:
                break;
            }
        }

    protected:
        JacobiPolynomialTpl<typename PS::NumScalar, PS::NOrder, PS::Options> jacobi_polynomial_;

    }; // class ControlParamModelJacobiPolynomialTpl

} // namespace galileo

#endif // __galileo_core_controls_control_param_jpoly_hpp__