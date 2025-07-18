#ifndef __galileo_multibody_contacts_contact_manager_hpp__
#define __galileo_multibody_contacts_contact_manager_hpp__

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "galileo/multibody/contacts/fwd.hpp"

#include "galileo/multibody/contacts/contact-base.hpp"
#include "galileo/multibody/force-base.hpp"

#include "galileo/multibody/contacts/contact-generic.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct ContactItemTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ContactItemTpl() {}
        ContactItemTpl(const std::string &name_, const Model_t &model_, const bool active_ = true)
            : name(name_), model(model_), active(active_) {}

        std::string name;
        Model_t model;
        bool active;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactManagerTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Collection_t = ContactCollectionTpl<PS>;
        using ModelManager_t = ContactModelManagerTpl<PS, ContactCollectionTpl>;
        using DataManager_t = ContactDataManagerTpl<PS, ContactCollectionTpl>;

        using Meta_t = ContactTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Item_t = ContactItemTpl<PS, ContactCollectionTpl>;

        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;

        using Jc_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNV_t::Value, PS::Options>;
        using a0_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1, PS::Options>;
        using da0_dx_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNDX_t::Value, PS::Options>;
        using dv_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, 1, PS::Options>;
        using ddv_dx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, PS::DimNDX_t::Value, PS::Options>;

        using Force_t = typename PS::Force_t;
        using ForceVector_t = GALILEO_ALIGNED_STD_VECTOR(Force_t);
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactDataManagerTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<ContactModelManagerTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    class ContactDataManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using Jc_t = typename traits<MetaManager_t>::Jc_t;
        using a0_t = typename traits<MetaManager_t>::a0_t;
        using da0_dx_t = typename traits<MetaManager_t>::da0_dx_t;
        using dv_t = typename traits<MetaManager_t>::dv_t;
        using ddv_dx_t = typename traits<MetaManager_t>::ddv_dx_t;
        using Force_t = typename traits<MetaManager_t>::Force_t;
        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;

        template <typename DataCollector>
        ContactDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : Jc(model_manager.get_nc_total(), model_manager.get_ps().get_nv()),
              a0(model_manager.get_nc_total()),
              da0_dx(model_manager.get_nc_total(), model_manager.get_ps().get_ndx()),
              dv(model_manager.get_ps().get_nv()),
              ddv_dx(model_manager.get_ps().get_nv(), model_manager.get_ps().get_ndx()),
              fext(model_manager.get_state()->get_robot().njoints, Force_t::Zero())
        {
            Jc.setZero();
            a0.setZero();
            da0_dx.setZero();
            dv.setZero();
            ddv_dx.setZero();
            for (typename ModelContainer_t::const_iterator
                     it = model_manager.getContacts().begin();
                 it != model_manager.getContacts().end(); ++it)
            {
                const Item_t &item = it->second;
                contacts.insert(
                    std::make_pair(item.name, item.model.createData(collector)));
            }
        }

        DataContainer_t contacts;
        ForceVector_t fext;

        Jc_t Jc;
        a0_t a0;
        da0_dx_t da0_dx;
        dv_t dv;
        ddv_dx_t ddv_dx;

    }; // class ContactDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    class ContactModelManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using MetaManager_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;
        using ForceIterator_t = typename ForceVector_t::iterator;

        ContactModelManagerTpl(const PS &ps, const std::shared_ptr<State_t> &state)
            : ps_(ps), state_(state),
              nc_active_dim_(DimensionTpl<Eigen::Dynamic>(0)), nc_total_dim_(DimensionTpl<Eigen::Dynamic>(0)),
              compute_all_contacts_(true)
        {
        }

        void addContact(const std::string &name, const Model_t &model, bool active = true)
        {
            std::pair<typename ModelContainer_t::iterator, bool> ret =
                contacts_.insert(std::make_pair(
                    name, Item_t(name, model, active)));
            if (ret.second == false)
            {
                std::cout << "Warning: we couldn't add the " << name
                          << " contact item, it already existed." << std::endl;
            }
            else if (active)
            {
                nc_active_dim_ += model.get_nc();
                nc_total_dim_ += model.get_nc();
                active_set_.insert(name);
            }
            else if (!active)
            {
                nc_total_dim_ += model.get_nc();
                inactive_set_.insert(name);
            }
        }

        void removeContact(const std::string &name)
        {
            typename ModelContainer_t::iterator it = contacts_.find(name);
            if (it != contacts_.end())
            {
                nc_active_dim_ -= it->second.model.get_nc();
                nc_total_dim_ -= it->second.model.get_nc();
                contacts_.erase(it);
                inactive_set_.erase(name);
            }
            else
            {
                std::cout << "Warning: we couldn't remove the " << name
                          << " contact item, it doesn't exist." << std::endl;
            }
        }

        void changeContactStatus(const std::string &name, bool active)
        {
            typename ModelContainer_t::iterator it = contacts_.find(name);
            if (it != contacts_.end())
            {
                if (active && !it->second.active)
                {
                    nc_active_dim_ += it->second.model.get_nc();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    nc_active_dim_ -= it->second.model.get_nc();
                    active_set_.erase(name);
                    inactive_set_.insert(name);
                    it->second.active = active;
                }
                it->second.active = active;
            }
            else
            {
                std::cout << "Warning: we couldn't change the status of the " << name
                          << " contact item, it doesn't exist." << std::endl;
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    auto nc_dim_i = m_i.model.get_nc_dim();
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calc(d_i, x.derived());
                        segment(data.a0, nc_accum_i, nc_dim_i) = d_i.a0();
                        block(data.Jc, nc_accum_i, 0, nc_dim_i, ps_.nv_dim) = d_i.Jc();
                    }
                    else
                    {
                        segment(data.a0, nc_accum_i, nc_dim_i).setZero();
                        block(data.Jc, nc_accum_i, 0, nc_dim_i, ps_.nv_dim).setZero();
                    }
                    nc_accum_i += nc_dim_i;
                }
            }
            else
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calc(d_i, x.derived());
                        auto nc_dim_i = m_i.model.get_nc_dim();
                        segment(data.a0, nc_accum_i, nc_dim_i) = d_i.a0();
                        block(data.Jc, nc_accum_i, 0, nc_dim_i, ps_.nv_dim) = d_i.Jc();
                        nc_accum_i += nc_dim_i;
                    }
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    auto nc_dim_i = m_i.model.get_nc_dim();
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calcDiff(d_i, x.derived());
                        block(data.da0_dx, nc_accum_i, 0, nc_dim_i, ps_.ndx_dim) = d_i.da0_dx();
                    }
                    else
                    {
                        block(data.da0_dx, nc_accum_i, 0, nc_dim_i, ps_.ndx_dim).setZero();
                    }
                    nc_accum_i += nc_dim_i;
                }
            }
            else
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calcDiff(d_i, x.derived());
                        auto nc_dim_i = m_i.model.get_nc_dim();
                        block(data.da0_dx, nc_accum_i, 0, nc_dim_i, ps_.ndx_dim) = d_i.da0_dx();
                        nc_accum_i += nc_dim_i;
                    }
                }
            }
        }

        template <typename VectorNvType>
        void updateAcceleration(DataManager_t &data, const Eigen::MatrixBase<VectorNvType> &dv) const
        {
            data.dv = dv.derived();
        }

        template <typename ForceVectorType>
        void updateForce(DataManager_t &data, const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            for (ForceIterator_t it = data.fext.begin(); it != data.fext.end(); ++it)
            {
                *it = PS::Force_t::Zero();
            }

            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    Data_t &d_i = it_d->second;
                    auto nc_dim_i = m_i.model.get_nc_dim();
                    if (m_i.active)
                    {
                        const Eigen::VectorBlock<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1>, Eigen::Dynamic> force_i =
                            segment(force, nc_accum_i, nc_dim_i);
                        m_i.model.updateForce(d_i, force_i);
                        const pinocchio::JointIndex joint =
                            state_->get_robot().frames[d_i.frame()].parent;
                        data.fext[joint] = d_i.fext();
                    }
                    else
                    {
                        m_i.model.setZeroForce(d_i);
                    }
                    nc_accum_i += nc_dim_i;
                }
            }
            else
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    Data_t &d_i = it_d->second;
                    if (m_i.active)
                    {
                        auto nc_dim_i = m_i.model.get_nc_dim();
                        const Eigen::VectorBlock<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1>, Eigen::Dynamic> force_i =
                            segment(force, nc_accum_i, nc_dim_i);
                        m_i.model.updateForce(d_i, force_i);
                        const pinocchio::JointIndex joint =
                            state_->get_robot().frames[d_i.frame()].parent;
                        data.fext[joint] = d_i.fext();
                        nc_accum_i += nc_dim_i;
                    }
                    else
                    {
                        m_i.model.setZeroForce(d_i);
                    }
                }
            }
        }

        template <typename MatrixNvNdxType>
        void updateAccelerationDiff(DataManager_t &data, const Eigen::MatrixBase<MatrixNvNdxType> &ddv_dx) const
        {
            data.ddv_dx = ddv_dx.derived();
        }

        template <typename MatrixNcNdxType, typename MatrixNcNduType>
        void updateForceDiff(DataManager_t &data, const Eigen::MatrixBase<MatrixNcNdxType> &df_dx, const Eigen::MatrixBase<MatrixNcNduType> &df_du) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::const_iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    Data_t &d_i = it_d->second;
                    auto nc_dim_i = m_i.model.get_nc_dim();
                    if (m_i.active)
                    {
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNDX_t::Value>> df_dx_i =
                            block(df_dx, nc_accum_i, 0, nc_dim_i, ps_.ndx_dim);
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNU_t::Value>> df_du_i =
                            block(df_du, nc_accum_i, 0, nc_dim_i, ps_.nu_dim);
                        m_i.model.updateForceDiff(d_i, df_dx_i, df_du_i);
                    }
                    else
                    {
                        m_i.model.setZeroForceDiff(d_i);
                    }
                    nc_accum_i += nc_dim_i;
                }
            }
            else
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    Data_t &d_i = it_d->second;
                    if (m_i.active)
                    {
                        auto nc_dim_i = m_i.model.get_nc_dim();
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNDX_t::Value>> df_dx_i =
                            block(df_dx, nc_accum_i, 0, nc_dim_i, ps_.ndx_dim);
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNU_t::Value>> df_du_i =
                            block(df_du, nc_accum_i, 0, nc_dim_i, ps_.nu_dim);
                        m_i.model.updateForceDiff(d_i, df_dx_i, df_du_i);
                        nc_accum_i += nc_dim_i;
                    }
                    else
                    {
                        m_i.model.setZeroForceDiff(d_i);
                    }
                }
            }
        }

        void updateRneaDiff(DataManager_t &data, RobotData_t &robot_data) const
        {
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::const_iterator it_d, end_d;
            for (it_m = contacts_.begin(), end_m = contacts_.end(),
                it_d = data.contacts.begin(), end_d = data.contacts.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                const Data_t &d_i = it_d->second;
                if (m_i.active)
                {
                    switch (m_i.model.get_type())
                    {
                    case pinocchio::ReferenceFrame::LOCAL:
                        break;
                    case pinocchio::ReferenceFrame::WORLD:
                    case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                        robot_data.dtau_dq += d_i.dtau_dq();
                        break;
                    }
                }
            }
        }
        template <typename DataCollector>
        DataManager_t createData(DataCollector *const collector) const
        {
            return DataManager_t(*this, collector);
        }

        const PS &get_ps() const
        {
            return ps_;
        }

        const std::shared_ptr<State_t> &get_state() const
        {
            return state_;
        }

        const ModelContainer_t &get_contacts() const
        {
            return contacts_;
        }

        const DimensionTpl<Eigen::Dynamic> &get_nc_active_dim() const
        {
            return nc_active_dim_;
        }

        int get_nc_active() const
        {
            return nc_active_dim_.value();
        }

        const DimensionTpl<Eigen::Dynamic> &get_nc_total_dim() const
        {
            return nc_total_dim_;
        }

        int get_nc_total() const
        {
            return nc_total_dim_.value();
        }

    protected:
        const PS &ps_;
        std::shared_ptr<State_t> state_;
        ModelContainer_t contacts_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

        DimensionTpl<Eigen::Dynamic> nc_active_dim_;
        DimensionTpl<Eigen::Dynamic> nc_total_dim_;

        bool compute_all_contacts_;

    }; // class ContactModelManagerTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_manager_hpp__
