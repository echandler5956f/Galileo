#ifndef __galileo_core_costs_cost_visitor_base_hpp__
#define __galileo_core_costs_cost_visitor_base_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"
#include "galileo/core/costs/cost-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Family tag for costs
        struct CostFamily
        {
        };

        // Trait specialization for Cost family
        template <>
        struct UnaryVisitorFamilyTraits<CostFamily>
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = CostModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = CostDataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = CostModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = CostDataBase<DataType, PS>;
        };

        // Cost-specific unary visitor base
        template <typename CostVisitorDerived, typename ReturnType = void>
        using CostUnaryVisitorBase = UnaryVisitorBase<CostFamily, CostVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_costs_cost_visitor_base_hpp__
