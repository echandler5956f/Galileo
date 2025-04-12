#ifndef __galileo_core_constraints_constraint_basic_visitors_hpp__
#define __galileo_core_constraints_constraint_basic_visitors_hpp__

#include "galileo/core/constraints/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous constraints

    // Constraint model visitors

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void constraint_calc_zeroth_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType>
    inline void constraint_calc_zeroth_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void constraint_calc_first_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType>
    inline void constraint_calc_first_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename LowerBoundType,
              typename UpperBoundType>
    inline void constraint_update_bounds(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        const Eigen::MatrixBase<LowerBoundType> &lb,
        const Eigen::MatrixBase<UpperBoundType> &ub);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline const typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t &constraint_lb(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline const typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t &constraint_ub(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline int constraint_ng(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline int constraint_nh(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model);

    // Constraint data visitors

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t &constraint_H(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t &constraint_Hx(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t &constraint_Hu(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::G_t &constraint_G(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gx_t &constraint_Gx(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gu_t &constraint_Gu(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

} // namespace galileo

#endif // __galileo_core_constraints_constraint_basic_visitors_hpp__