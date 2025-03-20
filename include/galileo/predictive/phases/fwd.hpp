#ifndef __galileo_predictive_phases_fwd_hpp__
#define __galileo_predictive_phases_fwd_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename VarScalar, typename NumScalar, int Options>
        struct PhaseCollectionDefaultTpl;
        
        using PhaseCollectionDefault = PhaseCollectionDefaultTpl<double, double, 0>;

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename, typename, int> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
        struct PhaseModelTpl;

        template <
            typename VarScalar,
            typename NumScalar,
            int Options,
            template <typename, typename, int> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
        struct PhaseDataTpl;

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_fwd_hpp__