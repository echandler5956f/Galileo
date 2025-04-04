#ifndef __galileo_multibody_contacts_contact_basic_visitors_hpp__
#define __galileo_multibody_contacts_contact_basic_visitors_hpp__

#include "galileo/multibody/contacts/fwd.hpp"

namespace galileo
{

    namespace multibody
    {

        // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous contacts

        // Contact model visitors

        template <typename PhaseSpec,
                  typename StateVectorType>
        inline void contact_calc_zeroth_order(
            const ContactModelTpl<PhaseSpec> &contact_model,
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x);

        template <typename PhaseSpec,
                  typename StateVectorType>
        inline void contact_calc_first_order(
            const ContactModelTpl<PhaseSpec> &contact_model,
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<StateVectorType> &x);

        // Contact data visitors

        template <typename PhaseSpec,
                  typename ForceVectorType>
        inline void contact_update_force(
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<ForceVectorType> &force);

        template <typename PhaseSpec,
                  typename MatrixNcNdxType,
                  typename MatrixNcNuType>
        inline void contact_update_force_diff(
            ContactDataTpl<PhaseSpec> &contact_data,
            const Eigen::MatrixBase<MatrixNcNdxType> &dF_dx,
            const Eigen::MatrixBase<MatrixNcNuType> &dF_du);

        template <typename PhaseSpec>
        inline void contact_set_zero_force(
            ContactDataTpl<PhaseSpec> &contact_data);

        template <typename PhaseSpec>
        inline void contact_set_zero_force_diff(
            ContactDataTpl<PhaseSpec> &contact_data);

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_basic_visitors_hpp__