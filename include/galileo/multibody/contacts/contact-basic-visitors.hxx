#ifndef __galileo_multibody_contacts_contact_basic_visitors_hxx__
#define __galileo_multibody_contacts_contact_basic_visitors_hxx__

#include <vector>

#include <boost/fusion/container/generation/make_vector.hpp>
#include "galileo/multibody/contacts/contact-unary-visitor.hpp"

#include "galileo/multibody/contacts/contact-basic-visitors.hpp"

#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace multibody
    {

        template <typename StateVectorType>
        struct ContactCalcZerothOrderVisitor
            : fusion::ContactUnaryVisitorBase<ContactCalcZerothOrderVisitor<StateVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>>;

            template <typename ContactModel>
            static void algo(
                const multibody::ContactModelBase<ContactModel> &contact_model,
                typename multibody::ContactDataBase<typename ContactModel::ContactDataDerived> &contact_data,
                const Eigen::MatrixBase<StateVectorType> &x)
            {
                contact_model.calc(contact_data, x.derived());
            }
        };

        template <typename PhaseSpec,
                  typename StateVectorType>
        inline void contact_calc_zeroth_order(
            const ContactModelTpl<PhaseSpec> &contact_model,
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            typedef ContactCalcZerothOrderVisitor<StateVectorType> Algo;

            Algo::run(contact_model, contact_data, typename Algo::ArgsType(x.derived()));
        }

        template <typename StateVectorType>
        struct ContactCalcFirstOrderVisitor
            : fusion::ContactUnaryVisitorBase<ContactCalcFirstOrderVisitor<StateVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>>;

            template <typename ContactModel>
            static void algo(
                const multibody::ContactModelBase<ContactModel> &contact_model,
                typename multibody::ContactDataBase<typename ContactModel::ContactDataDerived> &contact_data,
                const Eigen::MatrixBase<StateVectorType> &x)
            {
                contact_model.calcDiff(contact_data, x.derived());
            }
        };

        template <typename PhaseSpec,
                  typename StateVectorType>
        inline void contact_calc_first_order(
            const ContactModelTpl<PhaseSpec> &contact_model,
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x)
        {
            typedef ContactCalcFirstOrderVisitor<StateVectorType> Algo;

            Algo::run(contact_model, contact_data, typename Algo::ArgsType(x.derived()));
        }

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_basic_visitors_hxx__