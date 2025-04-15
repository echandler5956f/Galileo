#ifndef __galileo_core_constraints_constraint_data_base_hpp__
#define __galileo_core_constraints_constraint_data_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"
#include "galileo/core/constraints/constraint-model-base.hpp"

// We use traits rather than PhaseSpec,
// because each constraint model has its own NH and NG
#define GALILEO_CONSTRAINT_DATA_TYPEDEF(Constraint) \
    using H_t = typename traits<Constraint>::H_t;   \
    using Hx_t = typename traits<Constraint>::Hx_t; \
    using Hu_t = typename traits<Constraint>::Hu_t; \
    using G_t = typename traits<Constraint>::G_t;   \
    using Gx_t = typename traits<Constraint>::Gx_t; \
    using Gu_t = typename traits<Constraint>::Gu_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ConstraintDataBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        FORWARD_ACCESSOR(H_t, H);
        FORWARD_ACCESSOR(Hx_t, Hx);
        FORWARD_ACCESSOR(Hu_t, Hu);
        FORWARD_ACCESSOR(G_t, G);
        FORWARD_ACCESSOR(Gx_t, Gx);
        FORWARD_ACCESSOR(Gu_t, Gu);

    protected:
        inline ConstraintDataBase()
        {
        }

        inline ConstraintDataBase(const ConstraintDataBase &clone)
        {
            *this = clone;
        }

        inline ConstraintDataBase &operator=(const ConstraintDataBase &clone)
        {
            return *this;
        }

    }; // struct ConstraintDataBase

} // namespace galileo

#endif // __galileo_core_constraints_constraint_data_base_hpp__
