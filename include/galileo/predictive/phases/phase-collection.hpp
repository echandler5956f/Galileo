#ifndef __galileo_predictive_phases_phase_collection_hpp__
#define __galileo_predictive_phases_phase_collection_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename BasicSpec>
    struct PhaseCollectionDefaultTpl
    {
    public:
        using BS = BasicSpec;

        using PhaseModelVariant_t = boost::variant<PhaseModelVoid>; // TODO: Add phase models
        using PhaseDataVariant_t = boost::variant<PhaseDataVoid>;   // TODO: Add phase data
    };

    template <typename BasicSpec>
    using PhaseModelVariantTpl = PhaseCollectionDefaultTpl<BasicSpec>::PhaseModelVariant_t;

    template <typename BasicSpec>
    using PhaseDataVariantTpl = PhaseCollectionDefaultTpl<BasicSpec>::PhaseDataVariant_t;

} // namespace galileo

#endif // __galileo_predictive_phases_phase_collection_hpp__
