#ifndef __galileo_predictive_jumps_fwd_hpp__
#define __galileo_predictive_jumps_fwd_hpp__

#include "galileo/predictive/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct JumpModelImpulseFwdDynTpl;
    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct JumpDataImpulseFwdDynTpl;

} // namespace galileo

#endif // __galileo_predictive_jumps_fwd_hpp__
