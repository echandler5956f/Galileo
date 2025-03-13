#ifndef __galileo_core_constraints_constraint_model_base_hpp__
#define __galileo_core_constraints_constraint_model_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

#define GALILEO_CONSTRAINT_BASIC_TYPEDEF(Constraint)                                    \
    using Scalar = typename traits<Constraint>::Scalar;                                 \
    using VarScalar = typename traits<Constraint>::VarScalar;                           \
    static constexpr int Options = traits<Constraint>::Options;                                \
    using ConstraintModelDerived = typename traits<Constraint>::ConstraintModelDerived; \
    using ConstraintDataDerived = typename traits<Constraint>::ConstraintDataDerived;

#define GALILEO_CONSTRAINT_CONSTANTS(Constraint)

#define GALILEO_CONSTRAINT_MODEL_TYPEDEF(Constraint) \
    using ResidualModel_t = typename traits<Constraint>::ResidualModel_t;

#define GALILEO_CONSTRAINT_DATA_TYPEDEF(Constraint)                     \
    using ResidualData_t = typename traits<Constraint>::ResidualData_t; \
    using H_t = typename traits<Node>::H_t;                             \
    using Hx_t = typename traits<Node>::Hx_t;                           \
    using Hu_t = typename traits<Node>::Hu_t;                           \
    using G_t = typename traits<Node>::G_t;                             \
    using Gx_t = typename traits<Node>::Gx_t;                           \
    using Gu_t = typename traits<Node>::Gu_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived>
        class ConstraintModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintDerived = typename traits<Derived>::ConstraintDerived;
            GALILEO_CONSTRAINT_BASIC_TYPEDEF(ConstraintDerived);
            GALILEO_CONSTRAINT_CONSTANTS(ConstraintDerived);
            GALILEO_CONSTRAINT_MODEL_TYPEDEF(ConstraintDerived);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

        protected:
            inline ConstraintModelBase()
            {
            }

            inline ConstraintModelBase(const ConstraintModelBase &clone)
            {
                *this = clone;
            }

            inline ConstraintModelBase &operator=(const ConstraintModelBase &clone)
            {
                return *this;
            }

            ResidualModel_t residual_;

        }; // class ConstraintModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_base_hpp__
