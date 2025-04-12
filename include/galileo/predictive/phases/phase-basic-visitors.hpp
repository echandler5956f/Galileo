#ifndef __galileo_predictive_phases_phase_basic_visitors_hpp__
#define __galileo_predictive_phases_phase_basic_visitors_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include "galileo/core/states/state-base.hpp"

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
        const std::size_t &maxiter,
        const typename PhaseSpec::NumScalar &tol);

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline const typename PhaseSpec::SegmentModel_t &phase_segment_model(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model);

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline typename PhaseSpec::NumScalar phase_period(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model);

    // Phase data visitors

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline typename PhaseSpec::SegmentDataVector_t &phase_segment_data_vector(
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data);

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hpp__