#ifndef __galileo_multibody_contacts_contact_model_base_hpp__
#define __galileo_multibody_contacts_contact_model_base_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename Derived, typename PhaseSpec>
        struct ContactModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactDerived = typename traits<Derived>::ContactDerived;
            using ContactModelDerived = typename traits<ContactDerived>::ContactModelDerived;
            using ContactDataDerived = typename traits<ContactDerived>::ContactDataDerived;

            template <typename StateVectorType>
            void calc(ContactDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                derived().calc(data, x);
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                derived().calcDiff(data, x);
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataDerived &data,
                             const Eigen::MatrixBase<ForceVectorType> &f)
            {
                derived().updateForce(data, f);
            }

            template <typename JacobianXType, typename JacobianUType>
            void updateForceDiff(ContactDataDerived &data,
                                 const Eigen::MatrixBase<JacobianXType> &df_dx,
                                 const Eigen::MatrixBase<JacobianUType> &df_du)
            {
                data.df_dx = df_dx;
                data.df_du = df_du;
            }

            void setZeroForce(ContactDataDerived &data) const
            {
                data.f.setZero();
                data.fext.setZero();
            }

            void setZeroForceDiff(ContactDataDerived &data) const
            {
                data.df_dx.setZero();
                data.df_du.setZero();
            }

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
