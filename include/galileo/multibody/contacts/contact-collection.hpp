#ifndef __galileo_multibody_contacts_contact_collection_hpp__
#define __galileo_multibody_contacts_contact_collection_hpp__

#include "galileo/multibody/contacts/fwd.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    namespace multibody
    {

        template <typename PhaseSpec>
        struct ContactCollectionDefaultTpl
        {
            using PS = PhaseSpec;

            using ContactModelVariant = boost::variant<ContactModelVoid>; // TODO: add contact models

            using ContactDataVariant = boost::variant<ContactDataVoid>; // TODO: add contact data

        }; // struct ContactCollectionDefaultTpl

        template <typename PhaseSpec>
        using ContactModelVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ContactModelVariant;

        template <typename PhaseSpec>
        using ContactDataVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ContactDataVariant;

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_collection_hpp__