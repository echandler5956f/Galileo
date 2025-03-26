#ifndef __galileo_core_constraints_constraint_collection_hpp__
#define __galileo_core_constraints_constraint_collection_hpp__

#include "galileo/core/constraints/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace core
    {

        template <typename _VarScalar, typename _NumScalar, int _Options>
        struct ConstraintCollectionDefaultTpl
        {
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            using Options = _Options;

            using ConstraintModelVariant = boost::variant<>; // TODO: add constraint models

            using ConstraintDataVariant = boost::variant<>; // TODO: add constraint data

        }; // struct ConstraintCollectionDefaultTpl

        template <typename VarScalar, typename NumScalar, int Options>
        using ConstraintModelVariantTpl = ConstraintCollectionDefaultTpl<VarScalar, NumScalar, Options>::ConstraintModelVariant;

        template <typename VarScalar, typename NumScalar, int Options>
        using ConstraintDataVariantTpl = ConstraintCollectionDefaultTpl<VarScalar, NumScalar, Options>::ConstraintDataVariant;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_collection_hpp__