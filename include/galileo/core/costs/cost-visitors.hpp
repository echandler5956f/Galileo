#ifndef __galileo_core_costs_cost_visitors_hpp__
#define __galileo_core_costs_cost_visitors_hpp__

#include "galileo/core/costs/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous costs

    // Cost model visitors

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void cost_calc_zeroth_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType>
    inline void cost_calc_zeroth_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void cost_calc_first_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType>
    inline void cost_calc_first_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename DataCollector>
    inline CostDataTpl<PhaseSpec, CostCollectionTpl> cost_create_data(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        DataCollector *const collector);

    // Cost data visitors

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::L_t &cost_L(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lx_t &cost_Lx(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lu_t &cost_Lu(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxx_t &cost_Lxx(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxu_t &cost_Lxu(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Luu_t &cost_Luu(
        const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data);

} // namespace galileo

#endif // __galileo_core_costs_cost_visitors_hpp__
