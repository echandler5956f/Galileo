#ifndef __galileo_core_constraints_constraint_model_base_hpp__
#define __galileo_core_constraints_constraint_model_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

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

        using BoundVector_t = typename traits<Meta_t>::BoundVector_t;

        using DimNH_t = typename traits<Meta_t>::DimNH_t;
        using DimNG_t = typename traits<Meta_t>::DimNG_t;

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

        template <typename LowerBoundType, typename UpperBoundType>
        void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                          const Eigen::MatrixBase<UpperBoundType> &ub)
        {
            this->derived().updateBounds(lb.derived(), ub.derived());
        }

        const BoundVector_t &get_lb() const
        {
            return this->derived().get_lb();
        }

        const BoundVector_t &get_ub() const
        {
            return this->derived().get_ub();
        }

        const PS &get_ps() const
        {
            return ps_;
        }

        const int get_nh() const
        {
            if constexpr (DimNH_t::IsFixed)
            {
                return DimNH_t::Value;
            }
            else
            {
                return nh_dim_.value();
            }
        }

        const DimNH_t &get_nh_dim() const
        {
            return nh_dim_;
        }

        const int get_ng() const
        {
            if constexpr (DimNG_t::IsFixed)
            {
                return DimNG_t::Value;
            }
            else
            {
                return ng_dim_.value();
            }
        }

        const DimNG_t &get_ng_dim() const
        {
            return ng_dim_;
        }

    protected:
        inline ConstraintModelBase(const PS &ps, const DimNH_t &nh_dim, const DimNG_t &ng_dim)
            : ps_(ps), nh_dim_(nh_dim), ng_dim_(ng_dim)
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
        DimNG_t ng_dim_;

    }; // class ConstraintModelBase

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_base_hpp__
