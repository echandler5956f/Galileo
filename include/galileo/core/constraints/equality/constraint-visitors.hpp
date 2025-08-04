#ifndef __galileo_core_constraints_equality_constraint_visitors_hpp__
#define __galileo_core_constraints_equality_constraint_visitors_hpp__

#include "galileo/core/constraints/equality/fwd.hpp"

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
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank);

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
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename DataCollector>
    inline ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> constraint_create_data(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        DataCollector *const collector);

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    inline int constraint_get_nh(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model);

    // Constraint data visitors

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t constraint_H(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t constraint_Hx(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl>
    inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t constraint_Hu(
        const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data);

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_visitors_hpp__
