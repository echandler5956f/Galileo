#ifndef __galileo_core_controls_control_param_model_base_hpp__
#define __galileo_core_controls_control_param_model_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ControlParamModelBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using NumScalar = typename PS::NumScalar;

        using DimNU_t = typename PS::DimNU_t;
        using DimNW_t = typename PS::DimNW_t;
        using DimNOrder_t = typename PS::DimNOrder_t;

        template <typename ControlParamVectorType>
        void calc(Data_t &data, const NumScalar t, const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calc(data, t, w);
        }

        template <typename ControlParamVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<ControlParamVectorType> &w) const
        {
            this->derived().calcDiff(data, w);
        }

        template <typename ControlVectorType>
        void params(Data_t &data, const NumScalar t, const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().params(data, t, u);
        }

        template <typename ControlBoundVectorType, typename ControlParamBoundVectorType>
        void convertBounds(const Eigen::MatrixBase<ControlBoundVectorType> &u_lb,
                           const Eigen::MatrixBase<ControlBoundVectorType> &u_ub,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_lb,
                           const Eigen::MatrixBase<ControlParamBoundVectorType> &w_ub) const
        {
            this->derived().convertBounds(u_lb, u_ub, w_lb, w_ub);
        }

        template <AssignmentOp op = SETTO, typename InputMatrixType, typename OutputMatrixType>
        void multiplyByJacobian(Data_t &data,
                                const Eigen::MatrixBase<InputMatrixType> &A,
                                Eigen::MatrixBase<OutputMatrixType> &out) const
        {
            this->derived().template multiplyByJacobian<op>(data, A, out);
        }

        template <AssignmentOp op = SETTO, typename InputMatrixType, typename OutputMatrixType>
        void multiplyJacobianTransposeBy(Data_t &data,
                                         const Eigen::MatrixBase<InputMatrixType> &A,
                                         Eigen::MatrixBase<OutputMatrixType> &out) const
        {
            this->derived().template multiplyJacobianTransposeBy<op>(data, A, out);
        }

        Data_t createData(MemoryArena &arena) const { return this->derived().createData(arena); }

        const PS &get_ps() const { return ps_; }

    protected:
        inline ControlParamModelBase(const PS &ps) : ps_(ps) {}
        inline ControlParamModelBase(const ControlParamModelBase &clone) : ps_(clone.ps_) {}
        inline ControlParamModelBase &operator=(const ControlParamModelBase &clone)
        {
            ps_ = clone.ps_;
            return *this;
        }

        PS ps_;

    }; // class ControlParamModelBase

} // namespace galileo

#endif // __galileo_core_controls_control_param_model_base_hpp__
