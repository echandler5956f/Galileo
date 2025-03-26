#ifndef __galileo_core_integrators_integrator_model_base_hpp__
#define __galileo_core_integrators_integrator_model_base_hpp__

#include "galileo/core/integrators/integrator-base.hpp"

#define GALILEO_INTEGRATOR_BASIC_TYPEDEF(Integrator)                                    \
    using Scalar = typename traits<Integrator>::Scalar;                                 \
    using VarScalar = typename traits<Integrator>::VarScalar;                           \
    static constexpr int Options = traits<Integrator>::Options;                         \
    using IntegratorModelDerived = typename traits<Integrator>::IntegratorModelDerived; \
    using IntegratorDataDerived = typename traits<Integrator>::IntegratorDataDerived;

#define GALILEO_INTEGRATOR_CONSTANTS(Integrator) \
    static constexpr IntegratorTypes IntegratorType = traits<Integrator>::IntegratorType;

#define GALILEO_ERK_INTEGRATOR_MODEL_TYPEDEF(ERKIntegrator) \


#define GALILEO_ERK_INTEGRATOR_DATA_TYPEDEF(ERKIntegrator)


namespace galileo
{
    namespace core
    {

        enum class IntegratorTypes : int
        {
            EXPLICIT_RK = 0,
            IMPLICIT_RK = 1,
            LIFTED_IMPLICIT_RK = 2
        }; // enum class IntegratorTypes

    } // namespace core

} // namespace galileo

#endif // __galileo_core_integrators_integrator_model_base_hpp__
