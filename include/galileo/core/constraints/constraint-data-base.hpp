#ifndef __galileo_core_constraints_constraint_data_base_hpp__
#define __galileo_core_constraints_constraint_data_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"
#include "galileo/core/constraints/constraint-model-base.hpp"

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        struct ConstraintDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = typename traits<Derived>::ConstraintDerived;
            GALILEO_CONSTRAINT_BASIC_TYPEDEF(ConstraintDerived);
            GALILEO_CONSTRAINT_CONSTANTS(ConstraintDerived);
            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

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

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_data_base_hpp__
