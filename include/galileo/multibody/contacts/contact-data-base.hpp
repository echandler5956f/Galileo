#ifndef __galileo_multibody_contacts_contact_data_base_hpp__
#define __galileo_multibody_contacts_contact_data_base_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"
#include "galileo/multibody/contacts/contact-model-base.hpp"

#include "galileo/multibody/force-base.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct ContactDataBase : ForceDataBase<ContactDataBase<Derived>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ContactDerived = typename traits<Derived>::ContactDerived;
            GALILEO_CONTACT_BASIC_TYPEDEF(ContactDerived);
            GALILEO_CONTACT_CONSTANTS(ContactDerived);
            GALILEO_CONTACT_DATA_TYPEDEF(ContactDerived);

            // Jc_t Jc;
            // Fx_t Fx;
            // Fu_t Fu;

        protected:
            inline ContactDataBase()
            {
            }

            inline ContactDataBase(const ContactDataBase &clone)
            {
                *this = clone;
            }

            inline ContactDataBase &operator=(const ContactDataBase &clone)
            {
                return *this;
            }

        }; // struct ContactDataBase

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_data_base_hpp__
