#ifndef __galileo_core_controls_control - param - jpoly_hpp__
#define __galileo_core_controls_control -param - jpoly_hpp__

#include "galileo/core/controls/controls-param-base.hpp"
#include "galileo/math/polynomial.hpp"
namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  int NU,
                  int NDeg>
        struct ControlParamJacobiPolynomialTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NU,
                  int _NDeg>
        struct traits<ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>>
        {
            using JacobiPolynomialDerived = JacobiPolynomialTpl<_NumScalar, _NDeg, _Options>;

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NU = _NU;
            static constexpr int NDeg = _NDeg;
            static constexpr int NW = NU * NDeg;

            using ControlParamDataDerived = ControlParamDataJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;
            using ControlParamModelDerived = ControlParamModelJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;

            using U_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            using W_t = Eigen::Matrix<VarScalar, NDeg, 1, Options>;
            using Uw_t = Eigen::Matrix<VarScalar, NU, NW, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NU,
                  int _NDeg>
        struct traits<ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>>
        {
            using ControlParamDerived = ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;
            using VarScalar = traits<ControlParamDerived>::VarScalar;
            using NumScalar = traits<ControlParamDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NU,
                  int _NDeg>
        struct traits<ControlParamModelJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>>
        {
            using ControlParamDerived = ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;
            using VarScalar = traits<ControlParamDerived>::VarScalar;
            using NumScalar = traits<ControlParamDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NU,
                  int _NDeg>
        struct ControlParamDataJacobiPolynomialTpl : public ControlParamDataBase<ControlParamDataJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ControlParamDerived = ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;
            GALILEO_CONTROL_PARAM_BASIC_TYPEDEF(ControlParamDerived);
            GALILEO_CONTROL_PARAM_CONSTANTS(ControlParamDerived);
            GALILEO_CONTROL_PARAM_DATA_TYPEDEF(ControlParamDerived);

        }; // class ControlParamDataJacobiPolynomialTpl

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NU,
                  int _NDeg>
        class ControlParamModelJacobiPolynomialTpl : public ControlParamModelBase<ControlParamModelJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ControlParamDerived = ControlParamJacobiPolynomialTpl<_VarScalar, _NumScalar, _Options, _NU, _NDeg>;
            GALILEO_CONTROL_PARAM_BASIC_TYPEDEF(ControlParamDerived);
            GALILEO_CONTROL_PARAM_CONSTANTS(ControlParamDerived);
            GALILEO_CONTROL_PARAM_MODEL_TYPEDEF(ControlParamDerived);

            template <typename ControlParamVectorType>
            void calc(ControlParamDataDerived &data, const NumScalar &t,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                jacobi_polynomial_.barycentricInterpolation(t, w.reshaped(NW, NDeg), data.U.derived());
            }

            template <typename ControlParamVectorType>
            void calcDiff(ControlParamDataDerived &data,
                          const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                jacobi_polynomial_.barycentricInterpolationDiff(t, w.reshaped(NW, NDeg), data.dU_dw.derived());
            }

            template <typename ControlVectorType>
            void params(ControlParamDataDerived &data, const NumScalar &t,
                        const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                for (std::size_t i = 0; i < NDeg; ++i)
                {
                    data.w.segment(i * NU, NU) = u;
                }
            }

            template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
            void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                               const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                               const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                               const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
            {
                for (std::size_t i = 0; i < NDeg; ++i)
                {
                    w_lb.segment(i * NU, NU) = u_lb;
                    w_ub.segment(i * NU, NU) = u_ub;
                }
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void multiplyByJacobian(
                ControlParamDataDerived &data,
                const Eigen::MatrixBase<InputMatrixType> &A,
                Eigen::MatrixBase<OutputMatrixType> &out,
                const AssignmentOp op = setto) const
            {
                switch (op)
                {
                case setto:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(0, i * NU, NW, NU) = NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                case addto:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(0, i * NU, NW, NU) += NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                case rmfrom:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(0, i * NU, NW, NU) -= NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                default:
                    break;
                }
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void multiplyJacobianTransposeBy(
                ControlParamDataDerived &data,
                const Eigen::MatrixBase<InputMatrixType> &A,
                Eigen::MatrixBase<OutputMatrixType> &out,
                const AssignmentOp op = setto) const
            {
                switch (op)
                {
                case setto:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(i * NU, 0, NU, NW) = NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                case addto:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(i * NU, 0, NU, NW) += NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                case rmfrom:
                    for (std::size_t i = 0; i < NDeg; ++i)
                    {
                        out.block(i * NU, 0, NU, NW) -= NumScalar(data.dU_dW(0, i * NU)) * A;
                    }
                    break;
                default:
                    break;
                }
            }

        protected:
            JacobiPolynomialDerived jacobi_polynomial_;

        }; // class ControlParamModelJacobiPolynomialTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_controls_control-param-jpoly_hpp__