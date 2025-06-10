#ifndef __galileo_core_costs_cost_data_base_hpp__
#define __galileo_core_costs_cost_data_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

// We use traits rather than PhaseSpec,
// because each cost model has its own residual
#define GALILEO_COST_DATA_TYPEDEF(Cost)         \
    using L_t = typename traits<Cost>::L_t;     \
    using Lx_t = typename traits<Cost>::Lx_t;   \
    using Lu_t = typename traits<Cost>::Lu_t;   \
    using Lxx_t = typename traits<Cost>::Lxx_t; \
    using Lxu_t = typename traits<Cost>::Lxu_t; \
    using Luu_t = typename traits<Cost>::Luu_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct CostDataBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_COST_DATA_TYPEDEF(Meta_t);

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

} // namespace galileo

#endif // __galileo_core_costs_cost_data_base_hpp__
