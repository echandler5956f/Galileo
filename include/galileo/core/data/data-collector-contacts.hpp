#ifndef __galileo_core_data_data_collector_contacts_hpp__
#define __galileo_core_data_data_collector_contacts_hpp__

#include "galileo/multibody/contacts/contact-manager.hpp"
#include "galileo/multibody/contacts/fwd.hpp"

namespace galileo
{

    // Contact data mixin
    template <typename Derived, typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactDataMixinTpl
    {
        using PS = PhaseSpec;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        ContactDataMixinTpl(ContactDataManager_t *data)
            : contacts(data) {}

        ContactDataManager_t *contacts;

    }; // struct ContactDataMixinTpl

} // namespace galileo

#endif // __galileo_core_data_data_collector_contacts_hpp__
