#ifndef __galileo_predictive_phases_phase_collection_hpp__
#define __galileo_predictive_phases_phase_collection_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phases.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace predictive
    {

        template<typename PhaseSpec>
        struct PhaseCollectionDefaultTpl
        {
        public:
            using PS = PhaseSpec;

            using PhaseModelVariant = boost::variant<>; // TODO: Add phase models

            using PhaseDataVariant = boost::variant<>;   // TODO: Add phase data
        };

        template<typename PhaseSpec>
        using PhaseModelVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::PhaseModelVariant;

        template<typename PhaseSpec>
        using PhaseDataVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::PhaseDataVariant;

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_collection_hpp__
