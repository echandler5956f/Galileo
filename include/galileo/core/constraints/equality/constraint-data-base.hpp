#ifndef __galileo_core_constraints_equality_constraint_data_base_hpp__
#define __galileo_core_constraints_equality_constraint_data_base_hpp__

#include "galileo/core/constraints/equality/constraint-base.hpp"

#define GALILEO_CONSTRAINT_DATA_TYPEDEF(Constraint) \
    using H_t = typename traits<Constraint>::H_t; \
    using Hx_t = typename traits<Constraint>::Hx_t; \
    using Hu_t = typename traits<Constraint>::Hu_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ConstraintDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(Meta_t);

        FORWARD_ACCESSOR(H_t, H);
        FORWARD_ACCESSOR(Hx_t, Hx);
        FORWARD_ACCESSOR(Hu_t, Hu);

    protected:
        inline ConstraintDataBase() {}
        inline ConstraintDataBase(const ConstraintDataBase &clone) { *this = clone; }
        inline ConstraintDataBase &operator=(const ConstraintDataBase &clone) { return *this; }

    }; // struct ConstraintDataBase

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_data_base_hpp__
