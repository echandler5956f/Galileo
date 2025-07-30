#ifndef __galileo_predictive_phases_phase_visitors_hpp__
#define __galileo_predictive_phases_phase_visitors_hpp__

#include "galileo/predictive/phases/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous phases

    // Phase model visitors

    template <typename BasicSpec,
              template <typename BS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_zeroth_order(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws);

    template <typename BasicSpec,
              template <typename BS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_first_order(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws);

    template <typename BasicSpec,
              template <typename BS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_quasi_static(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        Eigen::MatrixBase<ControlParamMatrixType> &ws,
        const int maxiter, const typename BasicSpec::NumScalar tol);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline PhaseDataTpl<BasicSpec, PhaseCollectionTpl> phase_create_data(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model);

    // Phase data visitors

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNext_t phase_XNext_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t phase_XNextx_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t phase_XNextw_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t phase_L_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t phase_Lx_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t phase_Lw_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t phase_Lxx_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t phase_Lxw_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t phase_Lww_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t phase_H_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t phase_Hx_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t phase_Hw_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t phase_G_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t phase_Gx_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t phase_Gw_at_i(
        const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const int i);

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hpp__
