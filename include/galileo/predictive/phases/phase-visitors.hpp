#ifndef __galileo_predictive_phases_phase_visitors_hpp__
#define __galileo_predictive_phases_phase_visitors_hpp__

#include "galileo/predictive/phases/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous phases

    // Phase model visitors

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_zeroth_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws);

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_first_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws);

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_quasi_static(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        Eigen::MatrixBase<ControlParamMatrixType> &ws,
        const int maxiter, const typename PhaseSpec::NumScalar tol);

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename DataCollector>
    inline PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> phase_create_data(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        DataCollector *const collector);

    template <typename PhaseSpec, template <typename> class PhaseCollectionTpl>
    inline PhaseSpec phase_get_ps(const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model);

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hpp__
