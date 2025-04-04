#ifndef __galileo_multibody_contacts_contact_data_base_hpp__
#define __galileo_multibody_contacts_contact_data_base_hpp__

#include "galileo/multibody/contacts/contact-base.hpp"
#include "galileo/multibody/contacts/contact-model-base.hpp"

#define GALILEO_CONTACT_DATA_TYPEDEF(Contact)                                \
    using RobotDataPointer_t = typename traits<Contact>::RobotDataPointer_t; \
    using Index_t = typename traits<Contact>::Index_t;                       \
    using ReferenceFrame_t = typename traits<Contact>::ReferenceFrame_t;     \
    using SE3_t = typename traits<Contact>::SE3_t;                           \
    using MatrixNcNv_t = typename traits<Contact>::MatrixNcNv_t;             \
    using Force_t = typename traits<Contact>::Force_t;                       \
    using MatrixNcNdx_t = typename traits<Contact>::MatrixNcNdx_t;           \
    using MatrixNcNu_t = typename traits<Contact>::MatrixNcNu_t;             \
    using ActionMatrix_t = typename traits<Contact>::ActionMatrix_t;         \
    using MatrixNv_t = typename traits<Contact>::MatrixNv_t;

#define GALILEO_CONTACT_DATA_BASE_DEFAULT_ACCESSOR(Contact) \
    FORWARD_ACCESSOR(robot_data_pointer)                    \
    FORWARD_ACCESSOR(frame)                                 \
    FORWARD_ACCESSOR(type)                                  \
    FORWARD_ACCESSOR(jMf)                                   \
    FORWARD_ACCESSOR(Jc)                                    \
    FORWARD_ACCESSOR(f)                                     \
    FORWARD_ACCESSOR(fext)                                  \
    FORWARD_ACCESSOR(df_dx)                                 \
    FORWARD_ACCESSOR(df_du)                                 \
    FORWARD_ACCESSOR(fXj)                                   \
    FORWARD_ACCESSOR(a0)                                    \
    FORWARD_ACCESSOR(da0_dx)                                \
    FORWARD_ACCESSOR(dtau_dq)

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct ContactDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = typename traits<Derived>::PS;
            GALILEO_CONTACT_DATA_TYPEDEF(Derived)

            const RobotDataPointer_t &robot_data() const
            {
                return derived().robot_data_pointer_accessor();
            }

            RobotDataPointer_t &robot_data()
            {
                return derived().robot_data_pointer_accessor();
            }

            const Index_t &frame() const
            {
                return derived().frame_accessor();
            }

            Index_t &frame()
            {
                return derived().frame_accessor();
            }

            const ReferenceFrame_t &type() const
            {
                return derived().type_accessor();
            }

            ReferenceFrame_t &type()
            {
                return derived().type_accessor();
            }

            const SE3_t &jMf() const
            {
                return derived().jMf_accessor();
            }

            SE3_t &jMf()
            {
                return derived().jMf_accessor();
            }

            const MatrixNcNv_t &Jc() const
            {
                return derived().Jc_accessor();
            }

            MatrixNcNv_t &Jc()
            {
                return derived().Jc_accessor();
            }

            const Force_t &f() const
            {
                return derived().f_accessor();
            }

            Force_t &f()
            {
                return derived().f_accessor();
            }

            const Force_t &fext() const
            {
                return derived().fext_accessor();
            }

            Force_t &fext()
            {
                return derived().fext_accessor();
            }

            const MatrixNcNdx_t &df_dx() const
            {
                return derived().df_dx_accessor();
            }

            MatrixNcNdx_t &df_dx()
            {
                return derived().df_dx_accessor();
            }

            const MatrixNcNu_t &df_du() const
            {
                return derived().df_du_accessor();
            }

            MatrixNcNu_t &df_du()
            {
                return derived().df_du_accessor();
            }

            const ActionMatrix_t &fXj() const
            {
                return derived().fXj_accessor();
            }

            ActionMatrix_t &fXj()
            {
                return derived().fXj_accessor();
            }

            const VectorNc_t &a0() const
            {
                return derived().a0_accessor();
            }

            VectorNc_t &a0()
            {
                return derived().a0_accessor();
            }

            const MatrixNcNdx_t &da0_dx() const
            {
                return derived().da0_dx_accessor();
            }

            MatrixNcNdx_t &da0_dx()
            {
                return derived().da0_dx_accessor();
            }

            const MatrixNv_t &dtau_dq() const
            {
                return derived().dtau_dq_accessor();
            }

            MatrixNv_t &dtau_dq()
            {
                return derived().dtau_dq_accessor();
            }

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
