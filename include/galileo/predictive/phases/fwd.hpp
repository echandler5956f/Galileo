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

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct PhaseModelMultibodyTpl;

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct PhaseDataMultibodyTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class PhaseCollectionTpl>
    struct PhaseModelTpl;

    template <
        typename PhaseSpec,
        template <typename PS> class PhaseCollectionTpl>
    struct PhaseDataTpl;

} // namespace galileo

#endif // __galileo_predictive_phases_fwd_hpp__
