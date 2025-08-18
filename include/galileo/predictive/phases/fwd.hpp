#ifndef __galileo_predictive_phases_fwd_hpp__
#define __galileo_predictive_phases_fwd_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    struct PhaseModelVoid
    {
    }; // struct PhaseModelVoid`

    struct PhaseDataVoid
    {
    }; // struct PhaseDataVoid

    template <typename PhaseSpec, template <typename> class JumpTpl>
    struct PhaseModelDefaultTpl;
    template <typename PhaseSpec, template <typename> class JumpTpl>
    struct PhaseDataDefaultTpl;

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseModelTpl;
    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseDataTpl;

} // namespace galileo

#endif // __galileo_predictive_phases_fwd_hpp__
