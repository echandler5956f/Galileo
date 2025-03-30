#ifndef __galileo_core_costs_cost_data_base_hpp__
#define __galileo_core_costs_cost_data_base_hpp__

#include "galileo/core/costs/cost-base.hpp"
#include "galileo/core/costs/cost-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        struct CostDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using CostDerived = typename traits<Derived>::CostDerived;
            GALILEO_COST_BASIC_TYPEDEF(CostDerived);
            GALILEO_COST_CONSTANTS(CostDerived);
            GALILEO_COST_DATA_TYPEDEF(CostDerived);

            ResidualData_t residual;
            ActivationData_t activation;
            L_t L;
            Lx_t Lx;
            Lu_t Lu;
            Lxx_t Lxx;
            Lxu_t Lxu;
            Luu_t Luu;

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
