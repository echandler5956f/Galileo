#ifndef __galileo_predictive_phases_phase_basic_visitors_hpp__
#define __galileo_predictive_phases_phase_basic_visitors_hpp__

#include "galileo/predictive/phases/fwd.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline void calc_segment_from_phase(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &index);

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline void calc_diff_segment_from_phase(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &index);

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hpp__