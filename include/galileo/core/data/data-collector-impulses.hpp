#ifndef __galileo_core_data_data_collector_impulses_hpp__
#define __galileo_core_data_data_collector_impulses_hpp__

#include "galileo/multibody/impulses/impulse-manager.hpp"
#include "galileo/multibody/impulses/fwd.hpp"

namespace galileo
{

    // Impulse data mixin
    template <typename Derived, typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseDataMixinTpl
    {
        using PS = PhaseSpec;

        using ImpulseManagerMeta_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using ImpulseModelManager_t = typename traits<ImpulseManagerMeta_t>::ModelManager_t;
        using ImpulseDataManager_t = typename traits<ImpulseManagerMeta_t>::DataManager_t;

        ImpulseDataMixinTpl(std::shared_ptr<ImpulseDataManager_t> data)
            : impulses(data) {}

        std::shared_ptr<ImpulseDataManager_t> impulses;

    }; // struct ImpulseDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_impulses_hpp__
