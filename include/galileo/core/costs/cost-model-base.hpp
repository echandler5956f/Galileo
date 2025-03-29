#ifndef __galileo_core_costs_cost_model_base_hpp__
#define __galileo_core_costs_cost_model_base_hpp__

#include "galileo/core/costs/cost-base.hpp"

#define GALILEO_COST_BASIC_TYPEDEF(Residual, Activation)                                \
    using Scalar = typename traits<Residual>::Scalar;                                   \
    using VarScalar = typename traits<Residual>::VarScalar;                             \
    static constexpr int Options = traits<Residual>::Options;                           \
    using ResidualModelDerived = typename traits<Residual>::ResidualModelDerived;       \
    using ResidualDataDerived = typename traits<Residual>::ResidualDataDerived;         \
    using ActivationModelDerived = typename traits<Activation>::ActivationModelDerived; \
    using ActivationDataDerived = typename traits<Activation>::ActivationDataDerived;

#define GALILEO_COST_CONSTANTS(Cost)

#define GALILEO_COST_MODEL_TYPEDEF(Cost)                            \
    using ResidualModel_t = typename traits<Cost>::ResidualModel_t; \
    using ActivationModel_t = typename traits<Cost>::ActivationModel_t;

#define GALILEO_COST_DATA_TYPEDEF(Cost)                               \
    using ResidualData_t = typename traits<Cost>::ResidualData_t;     \
    using ActivationData_t = typename traits<Cost>::ActivationData_t; \
    using L_t = typename traits<Cost>::L_t;                           \
    using Lx_t = typename traits<Cost>::Lx_t;                         \
    using Lu_t = typename traits<Cost>::Lu_t;                         \
    using Lxx_t = typename traits<Cost>::Lxx_t;                       \
    using Lxu_t = typename traits<Cost>::Lxu_t;                       \
    using Luu_t = typename traits<Cost>::Luu_t;

namespace galileo
{
    namespace core
    {

        template <class _ResidualDerived, class _ActivationDerived>
        class CostModel
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ResidualDerived = typename traits<_ResidualDerived>::ResidualDerived;
            using ActivationDerived = typename traits<_ActivationDerived>::ActivationDerived;
            GALILEO_COST_BASIC_TYPEDEF(ResidualDerived, ActivationDerived);
            GALILEO_COST_CONSTANTS(ResidualDerived, ActivationDerived);
            GALILEO_COST_MODEL_TYPEDEF(ResidualDerived, ActivationDerived);

            template <typename StateVectorType, typename ControlVectorType>
            void calc(CostDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(CostDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            void calcCostDiff(CostDataDerived &cdata,
                              CostDataDerived &rdata,
                              const ActivationDataDerived &adata,
                              const bool update_u) const
            {
                derived().calcCostDiff(cdata, rdata, adata, update_u);
            }

        protected:
            ResidualModel_t residual_;
            ActivationModel_t activation_;

        }; // class CostModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_model_base_hpp__
