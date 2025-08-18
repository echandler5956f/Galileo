#ifndef __galileo_core_activations_fwd_hpp__
#define __galileo_core_activations_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{
    template <typename PhaseSpec, template <typename> class ResidualTpl>
    class ActivationModelQuadraticTpl;
    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct ActivationDataQuadraticTpl;

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    class ActivationModelWeightedQuadraticTpl;
    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct ActivationDataWeightedQuadraticTpl;

} // namespace galileo

#endif // __galileo_core_activations_fwd_hpp__
