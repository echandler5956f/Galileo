#ifndef __galileo_core_phase_collections_hpp__
#define __galileo_core_phase_collections_hpp__

#include "galileo/core/phase/fwd.hpp"
// #include "galileo/core/phase/phases.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename _Scalar, int _Options>
    struct PhaseCollectionDefaultTpl
    {
        using Scalar = _Scalar;
        enum
        {
            Options = _Options
        };

        using PhaseModelVariant = boost::variant<PhaseModelVoid>;
        using PhaseDataVariant = boost::variant<PhaseDataVoid>;
    };

    using PhaseModelVariant = typename PhaseCollectionDefault::PhaseModelVariant;
    using PhaseDataVariant = typename PhaseCollectionDefault::PhaseDataVariant;

} // namespace galileo

#endif // __galileo_core_phase_collections_hpp__
