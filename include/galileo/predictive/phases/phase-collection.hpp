#ifndef __galileo_predictive_phases_phase_collection_hpp__
#define __galileo_predictive_phases_phase_collection_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phases.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace predictive
    {

        template<typename _VarScalar,
                 typename _NumScalar,
                 int _Options>
        struct PhaseCollectionDefaultTpl
        {
        public:
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            using Options = _Options;

            using PhaseModelVariant = boost::variant<>; // TODO: Add phase models
            using PhaseDataVariant = boost::variant<>;   // TODO: Add phase data
        };

        using PhaseModelVariant =  PhaseCollectionDefault::PhaseModelVariant;
        using PhaseDataVariant = PhaseCollectionDefault::PhaseDataVariant;

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_collection_hpp__
