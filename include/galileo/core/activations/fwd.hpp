#ifndef __galileo_core_activations_fwd_hpp__
#define __galileo_core_activations_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationBoundsTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelQuadraticBarrierTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataQuadraticBarrierTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelWeightedQuadraticBarrierTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataWeightedQuadraticBarrierTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelQuadraticTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataQuadraticTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelQuadFlatExpTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataQuadFlatExpTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelQuadFlatLogTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataQuadFlatLogTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelWeightedQuadTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataWeightedQuadTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelSmooth1NormTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataSmooth1NormTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModelSmooth2NormTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationDataSmooth2NormTpl;

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        class ActivationModel2NormBarrierTpl;
        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename V, typename N, int O> class ResidualModelTpl>
        struct ActivationData2NormBarrierTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_activations_fwd_hpp__