#ifndef __galileo_core_costs_cost_data_base_hpp__
#define __galileo_core_costs_cost_data_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

#define GALILEO_COST_DATA_TYPEDEF(Cost)         \
    using L_t = typename traits<Cost>::L_t;     \
    using Lx_t = typename traits<Cost>::Lx_t;   \
    using Lu_t = typename traits<Cost>::Lu_t;   \
    using Lxx_t = typename traits<Cost>::Lxx_t; \
    using Lxu_t = typename traits<Cost>::Lxu_t; \
    using Luu_t = typename traits<Cost>::Luu_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        struct CostDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostDerived = typename traits<Derived>::CostDerived;
            using CostDataDerived = typename traits<CostDerived>::CostDataDerived;
            using CostModelDerived = typename traits<CostDerived>::CostModelDerived;

            GALILEO_COST_DATA_TYPEDEF(CostDerived);

            FORWARD_ACCESSOR(L_t, L);
            FORWARD_ACCESSOR(Lx_t, Lx);
            FORWARD_ACCESSOR(Lu_t, Lu);
            FORWARD_ACCESSOR(Lxx_t, Lxx);
            FORWARD_ACCESSOR(Lxu_t, Lxu);
            FORWARD_ACCESSOR(Luu_t, Luu);

        protected:
            inline CostDataBase()
            {
            }

            inline CostDataBase(const CostDataBase &clone)
            {
                *this = clone;
            }

            inline CostDataBase &operator=(const CostDataBase &clone)
            {
                return *this;
            }

        }; // struct CostDataBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_data_base_hpp__
