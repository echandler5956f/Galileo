#ifndef __galileo_multibody_spatial_contacts_contact_collection_hpp__
#define __galileo_multibody_spatial_contacts_contact_collection_hpp__

#include "galileo/domains/multibody/spatial/contacts/fwd.hpp"
#include "galileo/domains/multibody/spatial/contacts/impl/contact-3d.hpp"
#include "galileo/domains/multibody/spatial/contacts/impl/contact-6d.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename PhaseSpec>
    struct ContactCollectionDefaultTpl
    {
        using PS = PhaseSpec;

        using ContactModelVariant_t = boost::variant<ContactModel3dTpl<PS>, ContactModel6dTpl<PS>>;
        using ContactDataVariant_t = boost::variant<ContactData3dTpl<PS>, ContactData6dTpl<PS>>;

    }; // struct ContactCollectionDefaultTpl

    template <typename PhaseSpec>
    using ContactModelVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ContactModelVariant_t;

    template <typename PhaseSpec>
    using ContactDataVariantTpl = ContactCollectionDefaultTpl<PhaseSpec>::ContactDataVariant_t;

} // namespace galileo

#endif // __galileo_multibody_spatial_contacts_contact_collection_hpp__
