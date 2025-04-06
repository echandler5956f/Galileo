#ifndef __galileo_core_activations_fwd_hpp__
#define __galileo_core_activations_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationBoundsTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelQuadraticBarrierTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataQuadraticBarrierTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelWeightedQuadraticBarrierTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataWeightedQuadraticBarrierTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelQuadraticTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataQuadraticTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelQuadFlatExpTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataQuadFlatExpTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelQuadFlatLogTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataQuadFlatLogTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelWeightedQuadTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataWeightedQuadTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelSmooth1NormTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataSmooth1NormTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModelSmooth2NormTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationDataSmooth2NormTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        class ActivationModel2NormBarrierTpl;
        template <typename PhaseSpec,
                  template <typename PS> class ResidualTpl>
        struct ActivationData2NormBarrierTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_fwd_hpp__