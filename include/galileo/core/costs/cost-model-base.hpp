#ifndef __galileo_core_costs_cost_model_base_hpp__
#define __galileo_core_costs_cost_model_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        class CostModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostDerived = typename traits<Derived>::CostDerived;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(CostDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(CostDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(CostDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(CostDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return derived().createData(collector);
            }

        protected:
            inline CostModelBase()
            {
            }

            inline CostModelBase(const CostModelBase &clone)
            {
                *this = clone;
            }

            inline CostModelBase &operator=(const CostModelBase &clone)
            {
                return *this;
            }

        }; // class CostModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_model_base_hpp__
