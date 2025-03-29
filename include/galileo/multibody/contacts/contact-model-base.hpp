#ifndef __galileo_multibody_contacts_contact_model_base_hpp__
#define __galileo_multibody_contacts_contact_model_base_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

#define GALILEO_CONTACT_BASIC_TYPEDEF(Contact)                                 \
    using VarScalar = typename traits<Contact>::VarScalar;                     \
    using NumScalar = typename traits<Contact>::NumScalar;                     \
    static constexpr int Options = traits<Contact>::Options;                   \
    using ContactModelDerived = typename traits<Contact>::ContactModelDerived; \
    using ContactDataDerived = typename traits<Contact>::ContactDataDerived;

#define GALILEO_CONTACT_CONSTANTS(Contact)           \
    static constexpr int NX = traits<Contact>::NX;   \
    static constexpr int NU = traits<Contact>::NU;   \
    static constexpr int NDX = traits<Contact>::NDX; \
    static constexpr int NQ = traits<Contact>::NQ;   \
    static constexpr int NV = traits<Contact>::NV;   \
    static constexpr int NC = traits<Contact>::NC;

#define GALILEO_CONTACT_MODEL_TYPEDEF(Contact)

#define GALILEO_CONTACT_DATA_TYPEDEF(Contact)    \
    using Jc_t = typename traits<Contact>::Jc_t; \
    using Fx_t = typename traits<Contact>::Fx_t; \
    using Fu_t = typename traits<Contact>::Fu_t;

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct ContactModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ContactDerived = typename traits<Derived>::ContactDerived;
            GALILEO_CONTACT_BASIC_TYPEDEF(ContactDerived);
            GALILEO_CONTACT_CONSTANTS(ContactDerived);
            GALILEO_CONTACT_MODEL_TYPEDEF(ContactDerived);

        protected:
            inline ContactModelBase()
            {
            }

            inline ContactModelBase(const ContactModelBase &clone)
            {
                *this = clone;
            }

            inline ContactModelBase &operator=(const ContactModelBase &clone)
            {
                return *this;
            }

        }; // struct ContactModelBase

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_model_base_hpp__
