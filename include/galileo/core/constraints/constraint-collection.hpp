#ifndef __galileo_core_constraints_constraint_collection_hpp__
#define __galileo_core_constraints_constraint_collection_hpp__

#include "galileo/core/constraints/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec>
        struct ConstraintCollectionDefaultTpl
        {
            using PS = PhaseSpec;

            using ConstraintModelVariant = boost::variant<>; // TODO: add constraint models

            using ConstraintDataVariant = boost::variant<>; // TODO: add constraint data

        }; // struct ConstraintCollectionDefaultTpl

        template <typename PhaseSpec>
        using ConstraintModelVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::ConstraintModelVariant;

        template <typename PhaseSpec>
        using ConstraintDataVariantTpl = ConstraintCollectionDefaultTpl<PhaseSpec>::ConstraintDataVariant;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_collection_hpp__