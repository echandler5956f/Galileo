#ifndef __galileo_core_constraints_equality_constraint_unary_visitor_hpp__
#define __galileo_core_constraints_equality_constraint_unary_visitor_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"
#include "galileo/core/constraints/equality/constraint-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Family tag for constraints
        struct ConstraintFamily {};

        // Trait specialization for Constraint family
        template <>
        struct UnaryVisitorFamilyTraits<ConstraintFamily>
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = ConstraintModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = ConstraintDataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = ConstraintModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = ConstraintDataBase<DataType, PS>;
        };

        // Constraint-specific unary visitor base
        template <typename ConstraintVisitorDerived, typename ReturnType = void>
        using ConstraintUnaryVisitorBase = UnaryVisitorBase<ConstraintFamily, ConstraintVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_unary_visitor_hpp__
