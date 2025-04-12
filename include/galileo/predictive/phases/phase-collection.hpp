#ifndef __galileo_predictive_phases_phase_collection_hpp__
#define __galileo_predictive_phases_phase_collection_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phases.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct PhaseCollectionDefaultTpl
    {
    public:
        using PS = PhaseSpec;

        using ModelVariant_t = boost::variant<PhaseModelVoid>; // TODO: Add phase models
        using DataVariant_t = boost::variant<PhaseDataVoid>; // TODO: Add phase data
    };

    template <typename PhaseSpec>
    using PhaseModelVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::ModelVariant_t;

    template <typename PhaseSpec>
    using PhaseDataVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::DataVariant_t;

} // namespace galileo

#endif // __galileo_predictive_phases_phase_collection_hpp__
