#ifndef __galileo_core_residual_model_base_hpp__
#define __galileo_core_residual_model_base_hpp__

#include "galileo/core/residual/residual-base.hpp"

#define GALILEO_RESIDUAL_MODEL_TYPEDEF_GENERIC(Residual, TYPENAME)                \
    typedef TYPENAME traits<Residual>::Scalar Scalar;                             \
    typedef TYPENAME traits<Residual>::ResidualModelDerived ResidualModelDerived; \
    typedef TYPENAME traits<Residual>::ResidualDataDerived ResidualDataDerived;   \
    typedef TYPENAME traits<Residual>::ResidualModel_t ResidualModel_t;

#define GALILEO_RESIDUAL_TYPEDEF_TEMPLATE(Residual) \
    GALILEO_RESIDUAL_MODEL_TYPEDEF_GENERIC(Residual, typename)

#define GALILEO_RESIDUAL_CAST_TYPE_SPECIALIZATION(ResidualModelTpl) \
    template <typename Scalar, typename NewScalar>                  \
    struct CastType<NewScalar, ResidualModelTpl<Scalar>>            \
    {                                                               \
        typedef ResidualModelTpl<NewScalar> type;                   \
    }

namespace galileo
{
    template <typename Derived>
    class ResidualModelBase : CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using ResidualDerived = typename traits<Derived>::ResidualDerived;
        GALILEO_RESIDUAL_TYPEDEF_TEMPLATE(ResidualDerived);

        template <typename StateVectorType, typename ControlVectorType>
        void calc(ResidualDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs,
                  const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calc(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calc(ResidualDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calc(data, xs.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(ResidualDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs,
                      const Eigen::MatrixBase<ControlVectorType> &us) const
        {
            derived().calcDiff(data, xs.derived(), us.derived());
        }

        template <typename StateVectorType>
        void calcDiff(ResidualDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

    protected:
        ResidualModel_t residual_;

    }; // class ResidualModelBase

} // namespace galileo

#endif // __galileo_core_residual_model_base_hpp__
