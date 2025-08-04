#ifndef __galileo_multibody_contacts_contact_base_hpp__
#define __galileo_multibody_contacts_contact_base_hpp__

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/force-base.hpp"

#include "galileo/predictive/phases/phase-spec.hpp"

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

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        // Forward propogating the traits of the derived class. Basically, instead of looking at ContactBaseTpl,
        // ForceDataBase will look at the traits of the derived class.
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
    struct ContactDataBase
        : public ForceDataBase<ContactDataBase<Derived, PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Base = ForceDataBase<ContactDataBase<Derived, PS>, PS>;

        // Now we can access the traits of the derived class
        GALILEO_CONTACT_DATA_TYPEDEF(Meta_t);

        // Accessors required by ContactDataBase
        FORWARD_ACCESSOR(typename PS::ActionMatrix_t, fXj);
        FORWARD_ACCESSOR(VectorNc_t, a0);
        FORWARD_ACCESSOR(MatrixNcNdx_t, da0_dx);
        FORWARD_ACCESSOR(typename PS::MatrixNv_t, dtau_dq);

        // We have to override the CRTP derived() method that ForceDataBase inherits from internal::CRTP so that
        // we can return the derived object of ContactDataBase rather than ContactDataBase itself

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

        // Fulfill the contract for ForceDataBase and forward calls to the next derived class.
        // This uses the overridden derived() method to find the grandchild (e.g., ContactDataTpl).

        auto &robot_accessor() { return this->derived().robot_accessor(); }
        const auto &robot_accessor() const { return this->derived().robot_accessor(); }

        auto &frame_accessor() { return this->derived().frame_accessor(); }
        const auto &frame_accessor() const { return this->derived().frame_accessor(); }

        auto &type_accessor() { return this->derived().type_accessor(); }
        const auto &type_accessor() const { return this->derived().type_accessor(); }

        auto &jMf_accessor() { return this->derived().jMf_accessor(); }
        const auto &jMf_accessor() const { return this->derived().jMf_accessor(); }

        auto &f_accessor() { return this->derived().f_accessor(); }
        const auto &f_accessor() const { return this->derived().f_accessor(); }

        auto &fext_accessor() { return this->derived().fext_accessor(); }
        const auto &fext_accessor() const { return this->derived().fext_accessor(); }

        auto &Jc_accessor() { return this->derived().Jc_accessor(); }
        const auto &Jc_accessor() const { return this->derived().Jc_accessor(); }

        auto &df_dx_accessor() { return this->derived().df_dx_accessor(); }
        const auto &df_dx_accessor() const { return this->derived().df_dx_accessor(); }

        auto &df_du_accessor() { return this->derived().df_du_accessor(); }
        const auto &df_du_accessor() const { return this->derived().df_du_accessor(); }

    protected:
        inline ContactDataBase()
        {
        }

        inline ContactDataBase(const ContactDataBase &clone)
            : Base(clone)
        {
            *this = clone;
        }

        inline ContactDataBase &operator=(const ContactDataBase &clone)
        {
            return *this;
        }

    }; // struct ContactDataBase

    template <typename Derived, typename PhaseSpec>
    struct ContactModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        // ContactModelBase is one level up in the hierarchy from ContactDataBase
        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using RobotModel_t = typename PS::RobotModel_t;
        using RobotData_t = typename PS::RobotData_t;
        using FrameIndex_t = typename PS::FrameIndex_t;
        using ReferenceFrame_t = typename PS::ReferenceFrame_t;

        using DimNC_t = typename traits<Meta_t>::DimNC_t;

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        template <typename ForceVectorType>
        void updateForce(Data_t &data,
                         const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            this->derived().updateForce(data, force);
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiff(Data_t &data,
                             const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                             const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            this->derived().updateForceDiffImpl(data, df_dx, df_du);
        }

        template <typename MatrixNcNdxType, typename MatrixNcNuType>
        void updateForceDiffImpl(Data_t &data,
                                 const Eigen::MatrixBase<MatrixNcNdxType> &df_dx,
                                 const Eigen::MatrixBase<MatrixNcNuType> &df_du) const
        {
            data.df_dx = df_dx;
            data.df_du = df_du;
        }

        void setZeroForce(Data_t &data) const
        {
            this->derived().setZeroForceImpl(data);
        }

        void setZeroForceImpl(Data_t &data) const
        {
            data.f.setZero();
            data.fext.setZero();
        }

        void setZeroForceDiff(Data_t &data) const
        {
            this->derived().setZeroForceDiffImpl(data);
        }

        void setZeroForceDiffImpl(Data_t &data) const
        {
            data.df_dx.setZero();
            data.df_du.setZero();
        }

        Data_t createData(RobotData_t *const robot) const
        {
            return this->derived().createData(robot);
        }

        FrameIndex_t get_id() const
        {
            return this->derived().get_id_impl();
        }

        FrameIndex_t get_id_impl() const
        {
            return id_;
        }

        void set_id(const FrameIndex_t &id)
        {
            this->derived().set_id_impl(id);
        }

        void set_id_impl(const FrameIndex_t &id)
        {
            id_ = id;
        }

        ReferenceFrame_t get_type() const
        {
            return this->derived().get_type_impl();
        }

        ReferenceFrame_t get_type_impl() const
        {
            return type_;
        }

        void set_type(const ReferenceFrame_t &type)
        {
            this->derived().set_type_impl(type);
        }

        void set_type_impl(const ReferenceFrame_t &type)
        {
            type_ = type;
        }

        const DimNC_t &get_nc_dim() const
        {
            return this->derived().get_nc_dim_impl();
        }

        const DimNC_t &get_nc_dim_impl() const
        {
            return nc_dim_;
        }

        int get_nc() const
        {
            return this->derived().get_nc_impl();
        }

        int get_nc_impl() const
        {
            return nc_dim_.value();
        }

    protected:
        inline ContactModelBase(const FrameIndex_t &id, const ReferenceFrame_t &type, const DimNC_t &nc_dim)
            : nc_dim_(nc_dim), id_(id), type_(type)
        {
        }

        inline ContactModelBase(const ContactModelBase &clone)
            : nc_dim_(clone.nc_dim_), id_(clone.id_), type_(clone.type_)
        {
        }

        inline ContactModelBase &operator=(const ContactModelBase &clone)
        {
            nc_dim_ = clone.nc_dim_;
            id_ = clone.id_;
            type_ = clone.type_;
            return *this;
        }

        DimNC_t nc_dim_;
        FrameIndex_t id_;
        ReferenceFrame_t type_;

    }; // struct ContactModelBase

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_base_hpp__
