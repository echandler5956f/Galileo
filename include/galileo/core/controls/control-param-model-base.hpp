#ifndef __galileo_core_controls_control_param_model_base_hpp__
#define __galileo_core_controls_control_param_model_base_hpp__

#include "galileo/core/controls/control-param-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ControlParamModelBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = typename PS::ControlParamMeta_t;
        using Model_t = typename PS::ControlParamModel_t;
        using Data_t = typename PS::ControlParamData_t;

        template <typename ControlParamVectorType>
        void calc(Data_t &data, const PS::NumScalar t,
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
        void params(Data_t &data, const PS::NumScalar t,
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

        Data_t createData() const
        {
            return this->derived().createData();
        }

        const PS &get_ps() const
        {
            return ps_;
        }

        /**
         * @brief Return the dimension of the control space
         */
        const int get_nu() const
        {
            return this->derived().get_nu_impl();
        }

        const int get_nu_impl() const
        {
            if constexpr (PS::DimNU_t::IsFixed)
            {
                return PS::DimNU_t::Value;
            }
            else
            {
                return ps_.NU_dim.value();
            }
        }

        const PS::DimNU_t &NUDim() const
        {
            return ps_.NU_dim;
        }

        /**
         * @brief Return the order of the control parameterization
         */
        const int get_norder() const
        {
            return this->derived().get_norder_impl();
        }

        const int get_norder_impl() const
        {
            if constexpr (PS::DimNOrder_t::IsFixed)
            {
                return PS::DimNOrder_t::Value;
            }
            else
            {
                return ps_.NOrder_dim.value();
            }
        }

        const PS::DimNOrder_t &NOrderDim() const
        {
            return ps_.NOrder_dim;
        }

        /**
         * @brief Return the dimension of the control parameter space
         */
        const int get_nw() const
        {
            return this->derived().get_nw_impl();
        }

        const int get_nw_impl() const
        {
            if constexpr (PS::DimNW_t::IsFixed)
            {
                return PS::DimNW_t::Value;
            }
            else
            {
                return ps_.NW_dim.value();
            }
        }

        const PS::DimNW_t &NWDim() const
        {
            return ps_.NW_dim;
        }

    protected:
        inline ControlParamModelBase(const PS &ps) : ps_(ps)
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

        const PS &ps_;

    }; // class ControlParamModelBase

} // namespace galileo

#endif // __galileo_core_controls_control_param_model_base_hpp__
