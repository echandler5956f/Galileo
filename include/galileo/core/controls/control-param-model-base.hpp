#ifndef __galileo_core_controls_control_param_model_base_hpp__
#define __galileo_core_controls_control_param_model_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ControlParamModelBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename PS::ControlParamMeta_t;
        using Model_t = typename PS::ControlParamModel_t;
        using Data_t = typename PS::ControlParamData_t;

        template <typename ControlParamVectorType>
        void calc(Data_t &data, const typename PS::NumScalar t,
                  const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calc(data, t, w.derived());
        }

        template <typename ControlParamVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calcDiff(data, w.derived());
        }

        template <typename ControlVectorType>
        void params(Data_t &data, const typename PS::NumScalar t,
                    const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().params(data, t, u.derived());
        }

        template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
        void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                           const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
        {
            this->derived().convertBounds(u_lb.derived(), u_ub.derived(), w_lb.derived(), w_ub.derived());
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyByJacobian(
            Data_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            this->derived().multiplyByJacobian(data, A.derived(), out.derived(), op);
        }

        template <typename InputMatrixType, typename OutputMatrixType>
        void multiplyJacobianTransposeBy(
            Data_t &data,
            const Eigen::MatrixBase<InputMatrixType> &A,
            Eigen::MatrixBase<OutputMatrixType> &out,
            const AssignmentOp op = setto) const
        {
            this->derived().multiplyJacobianTransposeBy(data, A.derived(), out.derived(), op);
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

} // namespace galileo

#endif // __galileo_core_controls_control_param_model_base_hpp__
