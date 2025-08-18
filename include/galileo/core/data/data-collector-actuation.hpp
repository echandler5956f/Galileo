#ifndef __galileo_core_data_data_collector_actuation_hpp__
#define __galileo_core_data_data_collector_actuation_hpp__

#include "galileo/core/actuations/actuation-data-base.hpp"

namespace galileo
{

    // Actuation data mixin
    template <typename Derived, typename PhaseSpec>
    struct ActuationDataMixinTpl
    {
        using PS = PhaseSpec;

        using ActuationData_t = typename PS::ActuationData_t;

        ActuationDataMixinTpl(std::shared_ptr<ActuationData_t> data) : actuation(data) {}

        std::shared_ptr<ActuationData_t> actuation;

    }; // struct ActuationDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_actuation_hpp__
