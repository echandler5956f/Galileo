#ifndef __galileo_core_costs_cost_collection_hpp__
#define __galileo_core_costs_cost_collection_hpp__

#include "galileo/core/costs/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec>
        struct CostCollectionDefaultTpl
        {
            using PhaseSpec = PhaseSpec;

            using CostModelVariant = boost::variant<>; // TODO: add cost models

            using CostDataVariant = boost::variant<>; // TODO: add cost data

        }; // struct CostCollectionDefaultTpl

        template <typename PhaseSpec>
        using CostModelVariantTpl = CostCollectionDefaultTpl<PhaseSpec>::CostModelVariant;

        template <typename PhaseSpec>
        using CostDataVariantTpl = CostCollectionDefaultTpl<PhaseSpec>::CostDataVariant;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_collection_hpp__