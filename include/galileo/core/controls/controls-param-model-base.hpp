#ifndef __galileo_core_controls_controls_param_model_base_hpp__
#define __galileo_core_controls_controls_param_model_base_hpp__

#include "galileo/core/controls/controls-param-base.hpp"

#define GALILEO_CONTROL_PARAM_BASIC_TYPEDEF(ControlParam)       \
    using VarScalar = typename traits<ControlParam>::VarScalar; \
    using NumScalar = typename traits<ControlParam>::NumScalar; \
    static constexpr int Options = traits<ControlParam>::Options;

#define GALILEO_CONTROL_PARAM_CONSTANTS(ControlParam)   \
    static constexpr int NU = traits<ControlParam>::NU; \
    static constexpr int NW = traits<ControlParam>::NW;

#define GALILEO_CONTROL_PARAM_MODEL_TYPEDEF(ControlParam)

#define GALILEO_CONTROL_PARAM_DATA_TYPEDEF(ControlParam) \
    using U_t = typename traits<ControlParam>::U_t;      \
    using W_t = typename traits<ControlParam>::W_t;      \
    using Uw_t = typename traits<ControlParam>::Uw_t;

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

            template <typename ControlParamVectorType>
            void calc(ControlParamDataDerived &data, const NumScalar t,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                derived().calc(data, t, w.derived());
            }

            template <typename ControlParamVectorType>
            void calcDiff(ControlParamDataDerived &data,
                          const Eigen::MatrixBase<ControlParamVectorType> &w) const
            {
                derived().calcDiff(data, w.derived());
            }

            template <typename ControlVectorType>
            void params(ControlParamDataDerived &data, const NumScalar t,
                        const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().params(data, t, u.derived());
            }

            template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
            void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                               const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                               const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                               const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
            {
                derived().convertBounds(u_lb.derived(), u_ub.derived(), w_lb.derived(), w_ub.derived());
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
