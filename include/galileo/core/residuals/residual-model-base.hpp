#ifndef __galileo_core_residuals_residual_model_base_hpp__
#define __galileo_core_residuals_residual_model_base_hpp__

#include "galileo/core/residuals/residual-base.hpp"

#define GALILEO_RESIDUAL_BASIC_TYPEDEF(Residual)                                  \
    using Scalar = typename traits<Residual>::Scalar;                             \
    using VarScalar = typename traits<Residual>::VarScalar;                       \
    static constexpr int Options = traits<Residual>::Options;                     \
    using ResidualModelDerived = typename traits<Residual>::ResidualModelDerived; \
    using ResidualDataDerived = typename traits<Residual>::ResidualDataDerived;

#define GALILEO_RESIDUAL_CONSTANTS(Residual)        \
    static constexpr int NX = traits<Residual>::NX; \
    static constexpr int NU = traits<Residual>::NU; \
    static constexpr int NR = traits<Residual>::NR;

#define GALILEO_RESIDUAL_MODEL_TYPEDEF(Residual)

#define GALILEO_RESIDUAL_DATA_TYPEDEF(Residual)           \
    using R_t = typename traits<Residual>::R_t;           \
    using Rx_t = typename traits<Residual>::Rx_t;         \
    using Ru_t = typename traits<Residual>::Ru_t;         \
    using Arr_Rx_t = typename traits<Residual>::Arr_Rx_t; \
    using Arr_Ru_t = typename traits<Residual>::Arr_Ru_t;

namespace galileo
{

    template <typename Derived>
    class ResidualModelBase : internal::CRTP<ResidualModelBase<Derived>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using ResidualDerived = typename traits<Derived>::ResidualDerived;
        GALILEO_RESIDUAL_BASIC_TYPEDEF(ResidualDerived);
        GALILEO_RESIDUAL_CONSTANTS(ResidualDerived);
        GALILEO_RESIDUAL_MODEL_TYPEDEF(ResidualDerived);

        template <typename StateVectorType, typename ControlVectorType>
        void calc(ResidualDataDerived &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x.derived(), u.derived());
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(ResidualDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x.derived(), u.derived());
        }

        void calcCostDiff(CostDataDerived &cdata,
                          ResidualDataDerived &rdata,
                          const ActivationDataDerived &adata,
                          const bool update_u) const
        {
            this->derived().calcCostDiff(cdata, rdata, adata, update_u);
        }

    protected:
        inline ResidualModelBase()
        {
        }

        inline ResidualModelBase(const ResidualModelBase &clone)
        {
            *this = clone;
        }

        inline ResidualModelBase &operator=(const ResidualModelBase &clone)
        {
            return *this;
        }

    }; // class ResidualModelBase

} // namespace galileo

#endif // __galileo_core_residuals_residual_model_base_hpp__
