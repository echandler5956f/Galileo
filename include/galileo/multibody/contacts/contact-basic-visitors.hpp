#ifndef __galileo_multibody_contacts_contact_basic_visitors_hpp__
#define __galileo_multibody_contacts_contact_basic_visitors_hpp__

#include "galileo/multibody/contacts/fwd.hpp"

namespace galileo
{

    // These are basically all the visitors that we need in order to homogenously iterate over a set of heterogeneous contacts

    // Contact model visitors

    template <typename PhaseSpec,
              template <typename> class ContactCollectionTpl,
              typename DataCollector>
    inline ContactDataTpl<PhaseSpec, ContactCollectionTpl> contact_create_data(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        DataCollector *const collector);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    inline void contact_calc_zeroth_order(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename StateVectorType>
    inline void contact_calc_first_order(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<StateVectorType> &x);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename ForceVectorType>
    inline void contact_update_force(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<ForceVectorType> &force);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl,
              typename MatrixNcNdxType,
              typename MatrixNcNuType>
    inline void contact_update_force_diff(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data,
        const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
        const Eigen::MatrixBase<MatrixNcNuType> &df_du);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_zero_force(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_zero_force_diff(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline int contact_nc(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t contact_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline void contact_set_id(
        const ContactModelTpl<PhaseSpec, ContactCollectionTpl> &contact_model,
        const typename traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>::Index_t &id);

    // Contact data visitors

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::RobotDataPointer_t contact_robot_data_pointer(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Index_t contact_frame(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ReferenceFrame_t contact_type(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::SE3_t contact_jMf(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNv_t contact_Jc(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t contact_f(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::Force_t contact_fext(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t contact_df_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNu_t contact_df_du(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::ActionMatrix_t contact_fXj(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::VectorNc_t contact_a0(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNcNdx_t contact_da0_dx(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    inline typename traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>::MatrixNv_t contact_dtau_dq(
        const ContactDataTpl<PhaseSpec, ContactCollectionTpl> &contact_data);

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_basic_visitors_hpp__