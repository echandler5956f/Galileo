#ifndef __galileo_core_constraint_model_base_hpp__
#define __galileo_core_constraint_model_base_hpp__

#include "galileo/core/constraint/constraint-base.hpp"

#define GALILEO_CONSTRAINT_MODEL_TYPEDEF_GENERIC(Constraint, TYPENAME)                  \
    typedef TYPENAME traits<Constraint>::Scalar Scalar;                                 \
    typedef TYPENAME traits<Constraint>::ConstraintModelDerived ConstraintModelDerived; \
    typedef TYPENAME traits<Constraint>::ConstraintDataDerived ConstraintDataDerived; \
    typedef TYPENAME traits<Constraint>::ResidualModel_t ResidualModel_t; \

#define GALILEO_CONSTRAINT_TYPEDEF_TEMPLATE(Constraint) \
    GALILEO_CONSTRAINT_MODEL_TYPEDEF_GENERIC(Constraint, typename)

#define GALILEO_CONSTRAINT_CAST_TYPE_SPECIALIZATION(ConstraintModelTpl) \
    template <typename Scalar, typename NewScalar>                      \
    struct CastType<NewScalar, ConstraintModelTpl<Scalar>>              \
    {                                                                   \
        typedef ConstraintModelTpl<NewScalar> type;                     \
    }

namespace galileo
{
    template <typename Derived>
    class ConstraintModelBase : CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using ConstraintDerived = typename traits<Derived>::ConstraintDerived;
        GALILEO_CONSTRAINT_TYPEDEF_TEMPLATE(ConstraintDerived);

        template <typename StateVectorType, typename ControlVectorType>
        void calc(ConstraintDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs,
                  const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calc(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calc(ConstraintDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calc(data, xs.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs,
                      const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calcDiff(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calcDiff(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

        protected:
            ResidualModel_t residual_;
            

    }; // class ConstraintModelBase

} // namespace galileo

#endif // __galileo_core_constraint_model_base_hpp__
