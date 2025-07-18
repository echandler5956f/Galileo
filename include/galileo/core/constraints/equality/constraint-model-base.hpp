#ifndef __galileo_core_constraints_equality_constraint_model_base_hpp__
#define __galileo_core_constraints_equality_constraint_model_base_hpp__

#include "galileo/core/constraints/equality/constraint-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    class ConstraintModelBase
        : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using DimNH_t = typename traits<Meta_t>::DimNH_t;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x.derived());
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        const PS &get_ps() const
        {
            return this->derived().get_ps_impl();
        }

        const PS &get_ps_impl() const
        {
            return ps_;
        }

        const DimNH_t &get_nh_dim() const
        {
            return this->derived().get_nh_dim_impl();
        }

        const DimNH_t &get_nh_dim_impl() const
        {
            return nh_dim_;
        }

        const int get_nh() const
        {
            return this->derived().get_nh_impl();
        }

        const int get_nh_impl() const
        {
            return nh_dim_.value();
        }

    protected:
        inline ConstraintModelBase(const PS &ps, const DimNH_t &nh_dim)
            : ps_(ps), nh_dim_(nh_dim)
        {
        }

        inline ConstraintModelBase(const ConstraintModelBase &clone)
        {
            *this = clone;
        }

        inline ConstraintModelBase &operator=(const ConstraintModelBase &clone)
        {
            return *this;
        }

        const PS &ps_;
        DimNH_t nh_dim_;

    }; // class ConstraintModelBase

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_model_base_hpp__
