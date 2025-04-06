#ifndef __galileo_multibody_contacts_contact_base_hpp__
#define __galileo_multibody_contacts_contact_base_hpp__

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/force-base.hpp"

#define GALILEO_CONTACT_DATA_TYPEDEF(Contact)                        \
    using ActionMatrix_t = typename traits<Contact>::ActionMatrix_t; \
    using VectorNc_t = typename traits<Contact>::VectorNc_t;         \
    using MatrixNv_t = typename traits<Contact>::MatrixNv_t;

namespace galileo
{
    namespace multibody
    {

        template <typename Derived, typename PhaseSpec>
        struct ContactBaseTpl;

        // We are basically forward propogating the responsibility of filling the traits
        // to the derived CRTP class, because ContactDataBaseTpl/ContactModelBaseTpl are CRTP base classes
        template <typename Derived, typename PhaseSpec>
        struct traits<ContactBaseTpl<Derived, PhaseSpec>>
        {
            using PS = PhaseSpec;
            using ContactBase = ContactDataBaseTpl<Derived, PS>;
            using ForceDerived = ContactBase;
            using ContactDerived = Derived;

            // Retrieve the traits of the derived class
            GALILEO_FORCE_DATA_TYPEDEF(ContactDerived);
            GALILEO_CONTACT_DATA_TYPEDEF(ContactDerived);
        };

        template <typename Derived, typename PhaseSpec>
        struct traits<ContactDataBaseTpl<Derived, PhaseSpec>>
        {
            using PS = PhaseSpec;
            using ContactBase = ContactBaseTpl<Derived, PS>;
            using ContactDerived = typename traits<ContactBase>::ContactDerived;
        };

        template <typename Derived, typename PhaseSpec>
        struct traits<ContactModelBaseTpl<Derived, PhaseSpec>>
        {
            using PS = PhaseSpec;
            using ContactBase = ContactBaseTpl<Derived, PS>;
            using ContactDerived = typename traits<ContactBase>::ContactDerived;
        };

        template <typename Derived, typename PhaseSpec>
        struct ContactDataBaseTpl : public ForceDataBase<ContactDataBaseTpl<Derived, PhaseSpec>, PhaseSpec>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;
            using ContactBase = typename traits<ContactBaseTpl<Derived, PS>>::ContactBase;
            using ForceDerived = typename traits<ContactBase>::ForceDerived;

            GALILEO_FORCE_DATA_TYPEDEF(ForceDerived);
            GALILEO_CONTACT_DATA_TYPEDEF(ContactBase);

            FORWARD_ACCESSOR(RobotDataPointer_t, robot_data_pointer);
            FORWARD_ACCESSOR(Index_t, frame);
            FORWARD_ACCESSOR(ReferenceFrame_t, type);
            FORWARD_ACCESSOR(SE3_t, jMf);
            FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
            FORWARD_ACCESSOR(Force_t, f);
            FORWARD_ACCESSOR(Force_t, fext);
            FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
            FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

            FORWARD_ACCESSOR(ActionMatrix_t, fXj);
            FORWARD_ACCESSOR(VectorNc_t, a0);
            FORWARD_ACCESSOR(MatrixNcNdx_t, da0_dx);
            FORWARD_ACCESSOR(MatrixNv_t, dtau_dq);

        protected:
            inline ContactDataBaseTpl()
            {
            }

            inline ContactDataBaseTpl(const ContactDataBaseTpl &clone)
            {
                *this = clone;
            }

            inline ContactDataBaseTpl &operator=(const ContactDataBaseTpl &clone)
            {
                return *this;
            }

        }; // struct ContactDataBaseTpl

        template <typename Derived, typename PhaseSpec>
        struct ContactModelBaseTpl : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            // ContactModelBase is one level up in the hierarchy from ContactDataBase
            using ContactDerived = typename traits<Derived>::ContactDerived;
            using ContactModelDerived = typename traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = typename traits<ContactDerived>::ContactDataDerived;

            template <typename StateVectorType>
            void calc(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataDerived &data,
                             const Eigen::MatrixBase<ForceVectorType> &f)
            {
                derived().updateForce(data, f.derived());
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

            int nc() const
            {
                return derived().nc_impl();
            }

            int nc_impl() const
            {
                return traits<Derived>::NC;
            }

            Index_t id() const
            {
                return derived().id_impl();
            }

            void set_id(const Index_t &id)
            {
                derived().set_id_impl(id);
            }

        protected:
            inline ContactModelBaseTpl()
            {
            }

            inline ContactModelBaseTpl(const ContactModelBaseTpl &clone)
            {
                *this = clone;
            }

            inline ContactModelBaseTpl &operator=(const ContactModelBaseTpl &clone)
            {
                return *this;
            }

        }; // struct ContactModelBaseTpl
    }
}

#endif // __galileo_multibody_contacts_contact_base_hpp__
