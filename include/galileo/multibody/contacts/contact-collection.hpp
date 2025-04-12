#ifndef __galileo_multibody_contacts_contact_collection_hpp__
#define __galileo_multibody_contacts_contact_collection_hpp__

#include "galileo/multibody/contacts/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ContactCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ModelVariant_t = boost::variant<ContactModelVoid>; // TODO: add contact models
        using DataVariant_t = boost::variant<ContactDataVoid>; // TODO: add contact data

    }; // struct ContactCollectionDefaultTpl

    template <typename PhaseSpec>
    using ContactModelVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ModelVariant_t;

    template <typename PhaseSpec>
    using ContactDataVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::DataVariant_t;

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_collection_hpp__