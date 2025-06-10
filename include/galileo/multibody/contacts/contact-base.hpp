#ifndef __galileo_multibody_contacts_contact_base_hpp__
#define __galileo_multibody_contacts_contact_base_hpp__

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/force-base.hpp"

#define GALILEO_CONTACT_DATA_TYPEDEF(Contact) \
    GALILEO_FORCE_DATA_TYPEDEF(Contact);      \
    using VectorNc_t = typename traits<Contact>::VectorNc_t;

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

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        // Forward propogating the traits of the derived class
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);
    };

    template <typename Derived, typename PhaseSpec>
    struct traits<ContactDataBase<Derived, PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename Derived, typename PhaseSpec>
    struct traits<ContactModelBase<Derived, PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename Derived, typename PhaseSpec>
    struct ContactDataBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        // Now we can access the traits of the derived class
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        // Accessors required by ForceDataBase
        FORWARD_ACCESSOR(RobotData_t *, robot);
        FORWARD_ACCESSOR(FrameIndex_t, frame);
        FORWARD_ACCESSOR(ReferenceFrame_t, type);
        FORWARD_ACCESSOR(SE3_t, jMf);
        FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
        FORWARD_ACCESSOR(Force_t, f);
        FORWARD_ACCESSOR(Force_t, fext);
        FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
        FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

        // Accessors required by ContactDataBase
        FORWARD_ACCESSOR(ActionMatrix_t, fXj);
        FORWARD_ACCESSOR(VectorNc_t, a0);
        FORWARD_ACCESSOR(MatrixNcNdx_t, da0_dx);
        FORWARD_ACCESSOR(MatrixNv_t, dtau_dq);

        /**We have to override the CRTP derived() method to return the
         derived object because ForceDataBase is multi-level CRTP**/

        /** Return reference to this as derived object */
        inline Derived &derived() & noexcept
        {
            return *static_cast<Derived *>(this);
        }
        /** Return reference to this as derived object */
        inline const Derived &derived() const & noexcept
        {
            return *static_cast<Derived const *>(this);
        }
        /** Return reference to this as derived object, when this is rvalue */
        inline Derived &&derived() && noexcept
        {
            return std::move(*static_cast<Derived *>(this));
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

    template <typename Derived, typename PhaseSpec>
    struct ContactModelBase : public internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        // ContactModelBase is one level up in the hierarchy from ContactDataBase
        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using RobotModel_t = typename traits<Meta_t>::RobotModel_t;
        using FrameIndex_t = typename traits<Meta_t>::FrameIndex_t;
        using ReferenceFrame_t = typename traits<Meta_t>::ReferenceFrame_t;

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return this->derived().createData(collector);
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x.derived());
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x.derived());
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            this->derived().updateForce(data, force.derived());
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiff(Data_t &data,
                             const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                             const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            this->derived().updateForceDiffImpl(data, df_dx.derived(), df_du.derived());
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiffImpl(Data_t &data,
                                 const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                                 const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            data.df_dx() = df_dx;
            data.df_du() = df_du;
        }

        void setZeroForce(Data_t &data) const
        {
            this->derived().setZeroForceImpl(data);
        }

        void setZeroForceImpl(Data_t &data) const
        {
            data.f().setZero();
            data.fext().setZero();
        }

        void setZeroForceDiff(Data_t &data) const
        {
            this->derived().setZeroForceDiffImpl(data);
        }

        void setZeroForceDiffImpl(Data_t &data) const
        {
            data.df_dx().setZero();
            data.df_du().setZero();
        }

        const RobotModel_t *robot() const
        {
            return this->derived().robot_impl();
        }

        FrameIndex_t id() const
        {
            return this->derived().id_impl();
        }

        void set_id(const FrameIndex_t &id)
        {
            this->derived().set_id_impl(id);
        }

        ReferenceFrame_t type() const
        {
            return this->derived().type_impl();
        }

        void set_type(const ReferenceFrame_t &type)
        {
            this->derived().set_type_impl(type);
        }

        int nc() const
        {
            return this->derived().nc_impl();
        }

        int nu() const
        {
            return this->derived().nu_impl();
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
