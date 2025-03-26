#ifndef __galileo_core_activations_fwd_hpp__
#define __galileo_core_activations_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationBoundsTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelQuadraticBarrierTpl;
        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataQuadraticBarrierTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelWeightedQuadraticBarrierTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelQuadTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelQuadFlatExpTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataQuadFlatExpTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelQuadFlatLogTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataQuadFlatLogTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelWeightedQuadTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataWeightedQuadTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelSmooth1NormTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataSmooth1NormTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModelSmooth2NormTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationDataSmooth2NormTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        class ActivationModel2NormBarrierTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NR = -1>
        struct ActivationData2NormBarrierTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_fwd_hpp__