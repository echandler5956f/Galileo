#ifndef __galileo_core_constraints_constraint_model_base_hpp__
#define __galileo_core_constraints_constraint_model_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        class ConstraintModelBase : internal::CRTP<ConstraintModelBase<Derived, PhaseSpec>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintDerived = typename traits<Derived>::ConstraintDerived;
            using ConstraintDataDerived = typename traits<ConstraintDerived>::ConstraintDataDerived;
            using ConstraintModelDerived = typename traits<ConstraintDerived>::ConstraintModelDerived;

            using BoundVector_t = typename traits<ConstraintDerived>::BoundVector_t;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return derived().createData(collector);
            }

            template <typename LowerBoundType, typename UpperBoundType>
            void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                              const Eigen::MatrixBase<UpperBoundType> &ub)
            {
                derived().updateBounds(lb.derived(), ub.derived());
            }

            const BoundVector_t &lb() const
            {
                return derived().lb_impl();
            }

            const BoundVector_t &ub() const
            {
                return derived().ub_impl();
            }

            int ng() const
            {
                return derived().ng_impl();
            }

            int nh() const
            {
                return derived().nh_impl();
            }

            int ng_impl() const
            {
                return traits<ConstraintDerived>::NG;
            }

            int nh_impl() const
            {
                return traits<ConstraintDerived>::NH;
            }

        protected:
            inline ConstraintModelBase()
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

        }; // class ConstraintModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_base_hpp__
