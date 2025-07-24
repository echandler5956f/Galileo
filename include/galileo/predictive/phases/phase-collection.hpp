#ifndef __galileo_predictive_phases_phase_collection_hpp__
#define __galileo_predictive_phases_phase_collection_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct PhaseCollectionDefaultTpl
    {
    public:
        using PS = PhaseSpec;

        using PhaseModelVariant_t = boost::variant<PhaseModelVoid>; // TODO: Add phase models
        using PhaseDataVariant_t = boost::variant<PhaseDataVoid>;   // TODO: Add phase data
    };

    template <typename PhaseSpec>
    using PhaseModelVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::PhaseModelVariant_t;

    template <typename PhaseSpec>
    using PhaseDataVariantTpl = PhaseCollectionDefaultTpl<PhaseSpec>::PhaseDataVariant_t;

} // namespace galileo

#endif // __galileo_predictive_phases_phase_collection_hpp__
