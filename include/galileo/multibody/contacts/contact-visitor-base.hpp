#ifndef __galileo_multibody_contacts_contact_visitor_base_hpp__
#define __galileo_multibody_contacts_contact_visitor_base_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"
#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{
    namespace fusion
    {
        // Family tag for contacts
        struct ContactFamily
        {
        };

        // Trait specialization for Contact family
        template <>
        struct UnaryVisitorFamilyTraits<ContactFamily>
        {
            template <typename PS, template <typename> class CollectionTpl>
            using ModelTpl = ContactModelTpl<PS, CollectionTpl>;

            template <typename PS, template <typename> class CollectionTpl>
            using DataTpl = ContactDataTpl<PS, CollectionTpl>;

            template <typename ModelType, typename PS>
            using ModelBase = ContactModelBase<ModelType, PS>;

            template <typename DataType, typename PS>
            using DataBase = ContactDataBase<DataType, PS>;
        };

        // Contact-specific unary visitor base
        template <typename ContactVisitorDerived, typename ReturnType = void>
        using ContactUnaryVisitorBase = UnaryVisitorBase<ContactFamily, ContactVisitorDerived, ReturnType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_visitor_base_hpp__
