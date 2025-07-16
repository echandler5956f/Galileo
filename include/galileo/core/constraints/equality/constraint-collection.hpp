#ifndef __galileo_core_constraints_equality_constraint_collection_hpp__
#define __galileo_core_constraints_equality_constraint_collection_hpp__

#include "galileo/core/constraints/equality/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ConstraintCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ConstraintModelVariant_t = boost::variant<ConstraintModelVoid>; // TODO: add constraint models
        using ConstraintDataVariant_t = boost::variant<ConstraintDataVoid>;   // TODO: add constraint data

    }; // struct ConstraintCollectionDefaultTpl

    template <typename PhaseSpec>
    using ConstraintModelVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::ConstraintModelVariant_t;

    template <typename PhaseSpec>
    using ConstraintDataVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::ConstraintDataVariant_t;

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_collection_hpp__
