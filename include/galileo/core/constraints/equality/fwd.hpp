#ifndef __galileo_core_constraints_equality_fwd_hpp__
#define __galileo_core_constraints_equality_fwd_hpp__

#include "galileo/core/fwd.hpp"
#include <type_traits>

namespace galileo
{

    struct ConstraintModelVoid
    {
    }; // struct ConstraintModelVoid`

    struct ConstraintDataVoid
    {
    }; // struct ConstraintDataVoid

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct ConstraintModelResidualTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl>
    struct ConstraintDataResidualTpl;

    template <typename PhaseSpec>
    struct ConstraintCollectionDefaultTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintModelTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintDataTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintModelManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintDataManagerTpl;

} // namespace galileo

#endif // __galileo_core_constraints_equality_fwd_hpp__
