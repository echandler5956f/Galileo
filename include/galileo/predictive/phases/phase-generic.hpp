#ifndef __galileo_predictive_phases_phase_generic_hpp__
#define __galileo_predictive_phases_phase_generic_hpp__

#include "galileo/predictive/phases/phase-collection.hpp"
#include "galileo/predictive/phases/phase-basic-visitors.hxx"
#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename VarScalar, 
        typename NumScalar, 
        int Options, 
        template <typename V, typename N, int O> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
        struct PhaseTpl;
        using Phase = PhaseTpl<double, double, 0>;

        template <typename _VarScalar, typename _NumScalar, int _Options, 
        template <typename V, typename N, int O> class _PhaseCollectionTpl>
        struct traits<PhaseTpl<_VarScalar, _NumScalar, _Options, _PhaseCollectionTpl>>
        {
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;
            
            using PhaseCollectionTpl = _PhaseCollectionTpl;
        };


    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_generic_hpp__
