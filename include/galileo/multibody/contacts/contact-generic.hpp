#ifndef __galileo_multibody_contacts_contact_generic_hpp__
#define __galileo_multibody_contacts_contact_generic_hpp__

#include "galileo/multibody/contacts/contact-collection.hpp"
#include "galileo/multibody/contacts/contact-basic-visitors.hxx"

namespace galileo
{
    namespace multibody
    {

        template <typename PhaseSpec>
        struct ContactTpl;

        template <typename PhaseSpec>
        struct traits<ContactTpl<PhaseSpec>>
        {
            using PS = PhaseSpec;

            using ContactDataDerived = ContactDataTpl<PhaseSpec>;
            using ContactModelDerived = ContactModelTpl<PhaseSpec>;

            static constexpr int NC = Eigen::Dynamic;
        };

        template <typename PhaseSpec>
        struct traits<ContactDataTpl<PhaseSpec>>
        {
            using ContactDerived = ContactTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct traits<ContactModelTpl<PhaseSpec>>
        {
            using ContactDerived = ContactTpl<PhaseSpec>;
        };

        template <typename PhaseSpec>
        struct ContactDataTpl : public ContactDataBase<ContactDataTpl<PhaseSpec>>, PhaseSpec::ContactCollectionTpl<PhaseSpec>::ContactDataVariant
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
        };

        template <typename PhaseSpec>
        struct ContactModelTpl : public ContactModelBase<ContactModelTpl<PhaseSpec>>, PhaseSpec::ContactCollectionTpl<PhaseSpec>::ContactModelVariant
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactDerived = ContactTpl<PhaseSpec>;
            using ContactModelDerived = traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = traits<ContactDerived>::ContactDataDerived;

            template <typename StateVectorType>
            void calc(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                contact_calc_zeroth_order(*this, data, x);
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                contact_calc_first_order(*this, data, x);
            }

        }; // struct ContactModelTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_generic_hpp__
