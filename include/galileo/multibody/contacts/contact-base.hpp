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

    template <typename Derived, typename PhaseSpec>
    struct ContactBaseTpl;

    template <typename Derived, typename PhaseSpec>
    struct ContactDataBase;

    template <typename Derived, typename PhaseSpec>
    struct ContactModelBase;

    // We are basically forward propogating the responsibility of filling the traits
    // to the derived CRTP class, because ContactDataBase/ContactModelBase are CRTP base classes
    template <typename Derived, typename PhaseSpec>
    struct traits<ContactBaseTpl<Derived, PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactBaseTpl<Derived, PS>;
        using Data_t = ContactDataBase<Derived, PS>;
        using Model_t = ContactModelBase<Derived, PS>;

        // Retrieve the traits of the derived class
        GALILEO_FORCE_DATA_TYPEDEF(Meta_t);
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);
    };

    template <typename Derived, typename PhaseSpec>
    struct traits<ContactDataBase<Derived, PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactBaseTpl<Derived, PS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename Derived, typename PhaseSpec>
    struct traits<ContactModelBase<Derived, PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ContactBaseTpl<Derived, PS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename Derived, typename PhaseSpec>
    struct ContactDataBase : public ForceDataBase<ContactDataBase<Derived, PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ContactBaseTpl<Derived, PS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        GALILEO_FORCE_DATA_TYPEDEF(Meta_t);
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

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

    template <typename Derived, typename PhaseSpec>
    struct ContactModelBase : internal::CRTP<ContactModelBase<Derived, PhaseSpec>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        // ContactModelBase is one level up in the hierarchy from ContactDataBase
        using Meta_t = ContactBaseTpl<Derived, PS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using Index_t = typename traits<Meta_t>::Index_t;

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x)
        {
            this->derived().calc(data, x.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
        {
            this->derived().calcDiff(data, x.derived());
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force)
        {
            this->derived().updateForce(data, force.derived());
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiff(Data_t &data,
                             const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                             const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            this->derived().updateForceDiff(data, df_dx.derived(), df_du.derived());
        }

        void setZeroForce(Data_t &data) const
        {
            this->derived().setZeroForce(data);
        }

        void setZeroForceDiff(Data_t &data) const
        {
            this->derived().setZeroForceDiff(data);
        }

        int nc() const
        {
            return this->derived().nc_impl();
        }

        int nc_impl() const
        {
            return traits<Meta_t>::NC;
        }

        Index_t id() const
        {
            return this->derived().id_impl();
        }

        void set_id(const Index_t &id)
        {
            this->derived().set_id_impl(id);
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

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_base_hpp__
