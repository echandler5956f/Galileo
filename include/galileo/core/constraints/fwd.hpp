#ifndef __galileo_core_constraints_fwd_hpp__
#define __galileo_core_constraints_fwd_hpp__

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

    enum class ConstraintType
    {
        Equality = 0,
        Inequality = 1,
        Any = 2
    };

    template <ConstraintType EqualityInequality>
    inline constexpr bool is_equality_v = (EqualityInequality == ConstraintType::Equality);

    template <ConstraintType EqualityInequality>
    inline constexpr bool is_inequality_v = (EqualityInequality == ConstraintType::Inequality);

    template <ConstraintType EqualityInequality>
    inline constexpr bool is_any_constraint_v = (EqualityInequality == ConstraintType::Any);

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
    struct ConstraintModelResidualTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class ResidualTpl,
        ConstraintType EqualityInequality>
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

#endif // __galileo_core_constraints_fwd_hpp__
