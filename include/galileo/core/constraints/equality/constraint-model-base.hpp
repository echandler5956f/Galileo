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
            this->derived().calc(data, x, u);
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x, u);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        int get_nh() const
        {
            return this->derived().get_nh_impl();
        }

        int get_nh_impl() const
        {
            return nh_dim_.value();
        }

    protected:
        inline ConstraintModelBase(const DimNH_t &nh_dim)
            : nh_dim_(nh_dim)
        {
        }

        inline ConstraintModelBase(const ConstraintModelBase &clone)
            : nh_dim_(clone.nh_dim_)
        {
        }

        inline ConstraintModelBase &operator=(const ConstraintModelBase &clone)
        {
            nh_dim_ = clone.nh_dim_;
            return *this;
        }

        DimNH_t nh_dim_;

    }; // class ConstraintModelBase

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_model_base_hpp__
