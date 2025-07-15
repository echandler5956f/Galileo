#ifndef __galileo_core_costs_cost_collection_hpp__
#define __galileo_core_costs_cost_collection_hpp__

#include "galileo/core/costs/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct CostCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using CostModelVariant_t = boost::variant<CostModelVoid>; // TODO: add cost models
        using CostDataVariant_t = boost::variant<CostDataVoid>;   // TODO: add cost data

    }; // struct CostCollectionDefaultTpl

    template <typename PhaseSpec>
    using CostModelVariantTpl = CostCollectionDefaultTpl<PhaseSpec>::CostModelVariant_t;

    template <typename PhaseSpec>
    using CostDataVariantTpl = CostCollectionDefaultTpl<PhaseSpec>::CostDataVariant_t;

} // namespace galileo

#endif // __galileo_core_costs_cost_collection_hpp__
