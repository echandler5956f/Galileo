#ifndef __galileo_multibody_contacts_contact_collection_hpp__
#define __galileo_multibody_contacts_contact_collection_hpp__

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/contacts/implementations/contact-3d.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ContactCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ModelVariant_t = boost::variant<ContactModel3dTpl<PS>>; // TODO: add contact models
        using DataVariant_t = boost::variant<ContactData3dTpl<PS>>;   // TODO: add contact data

    }; // struct ContactCollectionDefaultTpl

    template <typename PhaseSpec>
    using ContactModelVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ModelVariant_t;

    template <typename PhaseSpec>
    using ContactDataVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::DataVariant_t;

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_collection_hpp__
