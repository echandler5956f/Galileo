#ifndef __galileo_core_constraints_fwd_hpp__
#define __galileo_core_constraints_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
    {

        struct ConstraintModelVoid
        {
        }; // struct ConstraintModelVoid`

        struct ConstraintDataVoid
        {
        }; // struct ConstraintDataVoid

        template <typename VarScalar, typename NumScalar, int Options>
        struct ConstraintCollectionDefaultTpl;

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct ConstraintModelTpl;

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename V, typename N, int O> class ConstraintCollectionTpl>
        struct ConstraintDataTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_fwd_hpp__