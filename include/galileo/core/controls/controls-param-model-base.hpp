#ifndef __galileo_core_controls_controls_param_model_base_hpp__
#define __galileo_core_controls_controls_param_model_base_hpp__

#include "galileo/core/controls/controls-param-base.hpp"

#define GALILEO_CONTROL_PARAM_BASIC_TYPEDEF(ControlParam)       \
    using NumScalar = typename traits<ControlParam>::NumScalar; \
    using VarScalar = typename traits<ControlParam>::VarScalar; \
    static constexpr int Options = traits<ControlParam>::Options;

#define GALILEO_CONTROL_PARAM_CONSTANTS(ControlParam)   \
    static constexpr int NU = traits<ControlParam>::NU; \
    static constexpr int NW = traits<ControlParam>::NW;

#define GALILEO_CONTROL_PARAM_MODEL_TYPEDEF(ControlParam)

#define GALILEO_CONTROL_PARAM_DATA_TYPEDEF(ControlParam) \
    using W_t = typename traits<ControlParam>::W_t;      \
    using U_t = typename traits<ControlParam>::U_t;      \
    using Wu_t = typename traits<ControlParam>::Wu_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class ControlParamModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ControlParamDerived = typename traits<Derived>::ControlParamDerived;
            GALILEO_CONTROL_PARAM_BASIC_TYPEDEF(ControlParamDerived);
            GALILEO_CONTROL_PARAM_CONSTANTS(ControlParamDerived);
            GALILEO_CONTROL_PARAM_MODEL_TYPEDEF(ControlParamDerived);

            template <typename ControlVectorType>
            void calc(ControlParamDataDerived &data, const NumScalar t,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, t, u.derived());
            }

            template <typename ControlVectorType>
            void calcDiff(ControlDataDerived &data,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, t, u.derived());
            }

            template <typename ControlParamVectorType>
            void params(ControlParamDataDerived &data, const NumScalar t,
                        const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                derived().params(data, t, w.derived());
            }

            template <typename ControlParamBoundVectorType, typename ControlBoundVectorType>
            void convertBounds(const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                               const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub,
                               const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                               const Eigen::MatrixBase<ControlBoundVectorType> &u_ub) const
            {
                derived().convertBounds(w_lb.derived(), w_ub.derived(), u_lb.derived(), u_ub.derived());
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void multiplyByJacobian(
                ControlParamDataDerived &data,
                const Eigen::MatrixBase<InputMatrixType> &A,
                Eigen::MatrixBase<OutputMatrixType> &out,
                const AssignmentOp op = setto) const
            {
                derived().multiplyByJacobian(data, A.derived(), out.derived(), op);
            }

            template <typename InputMatrixType, typename OutputMatrixType>
            void multiplyJacobianTransposeBy(
                ControlParamDataDerived &data,
                const Eigen::MatrixBase<InputMatrixType> &A,
                Eigen::MatrixBase<OutputMatrixType> &out,
                const AssignmentOp op = setto) const
            {
                derived().multiplyJacobianTransposeBy(data, A.derived(), out.derived(), op);
            }

        protected:
            inline ControlParamModelBase()
            {
            }

            inline ControlParamModelBase(const ControlParamModelBase &clone)
            {
                *this = clone;
            }

            inline ControlParamModelBase &operator=(const ControlParamModelBase &clone)
            {
                return *this;
            }

        }; // class ControlParamModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_controls_controls_param_model_base_hpp__
