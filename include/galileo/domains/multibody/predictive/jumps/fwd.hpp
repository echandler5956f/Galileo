#ifndef __galileo_multibody_predictive_jumps_fwd_hpp__
#define __galileo_multibody_predictive_jumps_fwd_hpp__

#include "galileo/domains/multibody/predictive/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct JumpModelImpulseFwdDynTpl;
    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct JumpDataImpulseFwdDynTpl;

} // namespace galileo

#endif // __galileo_multibody_predictive_jumps_fwd_hpp__
