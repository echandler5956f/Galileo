#ifndef __galileo_core_data_actuation_hpp__
#define __galileo_core_data_actuation_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/core/actuations/actuation-data-base.hpp"

namespace galileo
{

    namespace core
    {

        // Actuation data mixin
        template <typename Derived, typename PhaseSpec>
        struct ActuationDataMixinTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            using BS = typename PS::BS;

            std::shared_ptr<ActuationDataTpl<BS>> actuation;

            ActuationDataMixinTpl(std::shared_ptr<ActuationDataTpl<BS>> data)
                : actuation(data) {}

        }; // struct ActuationDataMixinTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_data_actuation_hpp__