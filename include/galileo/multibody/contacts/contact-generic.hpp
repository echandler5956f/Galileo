#ifndef __galileo_multibody_contacts_contact_generic_hpp__
#define __galileo_multibody_contacts_contact_generic_hpp__

#include "galileo/multibody/contacts/contact-collection.hpp"
#include "galileo/multibody/contacts/contact-basic-visitors.hxx"

namespace galileo
{
    namespace multibody
    {

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct ContactTpl;

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct traits<ContactTpl<PhaseSpec, ContactCollectionTpl>>
        {
            using PS = PhaseSpec;
            using ContactCollection = ContactCollectionTpl<PS>;

            static constexpr int NC = Eigen::Dynamic;

            using ContactDataDerived = ContactDataTpl<PS, ContactCollectionTpl>;
            using ContactModelDerived = ContactModelTpl<PS, ContactCollectionTpl>;

            // Traits required by ForceDataBase
            using RobotDataPointer_t = typename(PS::RobotData_t) *;
            using Index_t = typename PS::Index_t;
            using ReferenceFrame_t = typename PS::ReferenceFrame_t;
            using SE3_t = typename PS::SE3_t;
            using MatrixNcNv_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NV, PS::Options, 6, PS::NV>;
            using Force_t = typename PS::Force_t;
            using MatrixNcNdx_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NDX, PS::Options, 6, PS::NDX>;
            using MatrixNcNu_t = Eigen::Matrix<typename PS::VarScalar, NC, PS::NU, PS::Options, 6, PS::NU>;

            // Traits required by ContactDataBase
            using ActionMatrix_t = typename PS::ActionMatrix_t;
            using VectorNc_t = Eigen::Matrix<typename PS::VarScalar, NC, 1, PS::Options, 6, 1>;
            using MatrixNv_t = Eigen::Matrix<typename PS::VarScalar, PS::NV, PS::NV, PS::Options>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct traits<ContactDataTpl<PhaseSpec, ContactCollectionTpl>>
        {
            using PS = PhaseSpec;
            using ContactDerived = ContactTpl<PS, ContactCollectionTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct traits<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>
        {
            using PS = PhaseSpec;
            using ContactDerived = ContactTpl<PS, ContactCollectionTpl>;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct ContactDataTpl : public ContactDataBase<ContactDataTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>,
                                ContactCollectionTpl<PhaseSpec>::ContactDataVariant
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactDerived = ContactTpl<PS, ContactCollectionTpl>;
            using ContactModelDerived = typename traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = typename traits<ContactDerived>::ContactDataDerived;

            GALILEO_FORCE_DATA_TYPEDEF(ContactDerived);
            GALILEO_CONTACT_DATA_TYPEDEF(ContactDerived);

            using ContactCollection = ContactCollectionTpl<PS>;
            using ContactDataVariant = typename ContactCollection::ContactDataVariant;

            ContactDataVariant &toVariant()
            {
                return *static_cast<ContactDataVariant *>(this);
            }
            const ContactDataVariant &toVariant() const
            {
                return *static_cast<const ContactDataVariant *>(this);
            }

            RobotDataPointer_t robot_data_pointer() const
            {
                return galileo::multibody::contact_robot_data_pointer(*this);
            }

            Index_t frame() const
            {
                return galileo::multibody::contact_frame(*this);
            }

            ReferenceFrame_t type() const
            {
                return galileo::multibody::contact_type(*this);
            }

            SE3_t jMf() const
            {
                return galileo::multibody::contact_jMf(*this);
            }

            MatrixNcNv_t Jc() const
            {
                return galileo::multibody::contact_Jc(*this);
            }

            Force_t f() const
            {
                return galileo::multibody::contact_f(*this);
            }

            Force_t fext() const
            {
                return galileo::multibody::contact_fext(*this);
            }

            MatrixNcNdx_t df_dx() const
            {
                return galileo::multibody::contact_df_dx(*this);
            }

            MatrixNcNu_t df_du() const
            {
                return galileo::multibody::contact_df_du(*this);
            }

            ActionMatrix_t fXj() const
            {
                return galileo::multibody::contact_fXj(*this);
            }

            VectorNc_t a0() const
            {
                return galileo::multibody::contact_a0(*this);
            }

            MatrixNcNdx_t da0_dx() const
            {
                return galileo::multibody::contact_da0_dx(*this);
            }

            MatrixNv_t dtau_dq() const
            {
                return galileo::multibody::contact_dtau_dq(*this);
            }

            ContactDataTpl()
                : ContactDataVariant()
            {
            }

            ContactDataTpl(const ContactDataVariant &contact_data_variant)
                : ContactDataVariant(contact_data_variant)
            {
            }

            template <typename ConstraintDataDerived>
            ContactDataTpl(const ContactDataBase<ContactDataDerived> &contact_data)
                : ContactCollection::ContactDataVariant((ContactDataVariant)contact_data.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename ContactDataVariant::types, ContactDataDerived>));
            }

            GENERIC_ACCESSOR(RobotDataPointer_t, robot_data_pointer);
            GENERIC_ACCESSOR(Index_t, frame);
            GENERIC_ACCESSOR(ReferenceFrame_t, type);
            GENERIC_ACCESSOR(SE3_t, jMf);
            GENERIC_ACCESSOR(MatrixNcNv_t, Jc);
            GENERIC_ACCESSOR(Force_t, f);
            GENERIC_ACCESSOR(Force_t, fext);
            GENERIC_ACCESSOR(MatrixNcNdx_t, df_dx);
            GENERIC_ACCESSOR(MatrixNcNu_t, df_du);

            GENERIC_ACCESSOR(ActionMatrix_t, fXj);
            GENERIC_ACCESSOR(VectorNc_t, a0);
            GENERIC_ACCESSOR(MatrixNcNdx_t, da0_dx);
            GENERIC_ACCESSOR(MatrixNv_t, dtau_dq);
        };

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct ContactModelTpl : public ContactModelBase<ContactModelTpl<PhaseSpec, ContactCollectionTpl>>,
                                 ContactCollectionTpl<PhaseSpec>::ContactModelVariant
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactDerived = ContactTpl<PhaseSpec>;
            using ContactModelDerived = traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = traits<ContactDerived>::ContactDataDerived;

            using ContactCollection = ContactCollectionTpl<PS>;
            using ContactModelVariant = typename ContactCollection::ContactModelVariant;

            using Index_t = typename traits<ContactDerived>::Index_t;

            ContactModelTpl()
                : ContactModelVariant()
            {
            }

            ContactModelTpl(const ContactModelVariant &contact_model_variant)
                : ContactModelVariant(contact_model_variant)
            {
            }

            template <typename ContactModelDerived>
            ContactModelTpl(const ContactModelBase<ContactModelDerived> &contact_model)
                : ContactCollection::ContactModelVariant((ContactModelVariant)contact_model.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename ContactModelVariant::types, ContactModelDerived>));
            }

            ConstraintModelVariant &toVariant()
            {
                return *static_cast<ContactModelVariant *>(this);
            }

            const ConstraintModelVariant &toVariant() const
            {
                return *static_cast<const ContactModelVariant *>(this);
            }

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

            template <typename ForceVectorType>
            void updateForce(ContactDataDerived &data,
                             const Eigen::MatrixBase<ForceVectorType> &f)
            {
                contact_update_force(*this, data, f.derived());
            }

            template <typename JacobianXType, typename JacobianUType>
            void updateForceDiff(ContactDataDerived &data,
                                 const Eigen::MatrixBase<JacobianXType> &df_dx,
                                 const Eigen::MatrixBase<JacobianUType> &df_du)
            {
                data.df_dx() = df_dx;
                data.df_du() = df_du;
            }

            void setZeroForce(ContactDataDerived &data) const
            {
                data.f().setZero();
                data.fext().setZero();
            }

            void setZeroForceDiff(ContactDataDerived &data) const
            {
                data.df_dx().setZero();
                data.df_du().setZero();
            }

            int nc_impl() const
            {
                return galileo::multibody::contact_nc(*this);
            }

            Index_t id_impl() const
            {
                return galileo::multibody::contact_id(*this);
            }

            void set_id_impl(const Index_t &id)
            {
                galileo::multibody::contact_set_id(*this, id);
            }

        }; // struct ContactModelTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_generic_hpp__
