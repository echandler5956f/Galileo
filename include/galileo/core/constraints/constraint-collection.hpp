#ifndef __galileo_core_constraints_constraint_collection_hpp__
#define __galileo_core_constraints_constraint_collection_hpp__

#include "galileo/core/constraints/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ConstraintCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ModelVariant_t = boost::variant<ConstraintModelVoid>; // TODO: add constraint models
        using DataVariant_t = boost::variant<ConstraintDataVoid>;   // TODO: add constraint data

    }; // struct ConstraintCollectionDefaultTpl

    template <typename PhaseSpec>
    using ConstraintModelVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::ModelVariant_t;

    template <typename PhaseSpec>
    using ConstraintDataVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::DataVariant_t;

} // namespace galileo

#endif // __galileo_core_constraints_constraint_collection_hpp__
