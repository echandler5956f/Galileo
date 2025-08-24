#ifndef __galileo_core_data_data_collector_contacts_hpp__
#define __galileo_core_data_data_collector_contacts_hpp__

#include "galileo/domains/multibody/spatial/contacts/contact-manager.hpp"

namespace galileo
{

    // Contact data mixin
    template <typename Derived, typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct ContactDataMixinTpl
    {
        using PS = PhaseSpec;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        ContactDataMixinTpl(std::shared_ptr<ContactDataManager_t> data) : contacts(data) {}

        std::shared_ptr<ContactDataManager_t> contacts;

    }; // struct ContactDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_contacts_hpp__
