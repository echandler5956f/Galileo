#ifndef __galileo_multibody_contacts_contact_manager_hpp__
#define __galileo_multibody_contacts_contact_manager_hpp__

#include <iostream>
#include <string>
#include <map>
#include <set>

#include "galileo/multibody/contacts/fwd.hpp"

#include "galileo/multibody/force-base.hpp"
#include "galileo/multibody/contacts/contact-base.hpp"

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
        ContactItemTpl(const std::string &name_in, const Model_t &model_in, const bool active_in = true)
            : name(name_in), model(model_in), active(active_in) {}

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

        using Jc_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NV, PS::Options>;
        using a0_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1, PS::Options>;
        using da0_dx_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX, PS::Options>;
        using dv_t = Eigen::GMatrix<typename PS::VarScalar, PS::NV, 1, PS::Options>;
        using ddv_dx_t = Eigen::GMatrix<typename PS::VarScalar, PS::NV, PS::NDX, PS::Options>;

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

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

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
        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;

        Jc_t Jc;
        a0_t a0;
        da0_dx_t da0_dx;
        dv_t dv;
        ddv_dx_t ddv_dx;

        DataContainer_t contacts;
        ForceVector_t fext;

        template <typename DataCollector>
        ContactDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : Jc(model_manager.nc_total(), NV),
              a0(model_manager.nc_total()),
              da0_dx(model_manager.nc_total(), NDX),
              dv(NV),
              ddv_dx(NV, NDX),
              fext(model_manager.getState().get_robot()->njoints, Force_t::Zero())
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

    }; // class ContactDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    class ContactModelManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

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

        ContactModelManagerTpl(State_t *state) : state_(state), nc_(0), nc_total_(0), compute_all_contacts_(true) {}

        template <typename DataCollector>
        DataManager_t createData(DataCollector *const collector) const
        {
            return DataManager_t(*this, collector);
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
                nc_ += model.nc();
                nc_total_ += model.nc();
                active_set_.insert(name);
            }
            else if (!active)
            {
                nc_total_ += model.nc();
                inactive_set_.insert(name);
            }
        }

        void removeContact(const std::string &name)
        {
            typename ModelContainer_t::iterator it = contacts_.find(name);
            if (it != contacts_.end())
            {
                nc_ -= it->second.model.nc();
                nc_total_ -= it->second.model.nc();
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
                    nc_ += it->second.model.nc();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    nc_ -= it->second.model.nc();
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
            int nc = 0;
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    const int nc_i = m_i.model.nc();
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calc(d_i, x.derived());
                        data.a0.segment(nc, nc_i) = d_i.a0();
                        data.Jc.block(nc, 0, nc_i, PS::NV) = d_i.Jc();
                    }
                    else
                    {
                        data.a0.segment(nc, nc_i).setZero();
                        data.Jc.block(nc, 0, nc_i, PS::NV).setZero();
                    }
                    nc += nc_i;
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
                        const int nc_i = m_i.model.nc();
                        data.a0.segment(nc, nc_i) = d_i.a0();
                        data.Jc.block(nc, 0, nc_i, PS::NV) = d_i.Jc();
                        nc += nc_i;
                    }
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            int nc = 0;
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            if (compute_all_contacts_)
            {
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const Item_t &m_i = it_m->second;
                    const int nc_i = m_i.model.nc();
                    if (m_i.active)
                    {
                        Data_t &d_i = it_d->second;

                        m_i.model.calcDiff(d_i, x.derived());
                        data.da0_dx.block(nc, 0, nc_i, PS::NDX) = d_i.da0_dx();
                    }
                    else
                    {
                        data.da0_dx.block(nc, 0, nc_i, PS::NDX).setZero();
                    }
                    nc += nc_i;
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
                        const int nc_i = m_i.model.nc();
                        data.da0_dx.block(nc, 0, nc_i, PS::NDX) = d_i.da0_dx();
                        nc += nc_i;
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

            std::size_t nc = 0;
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
                    const int nc_i = m_i.model.nc();
                    if (m_i.active)
                    {
                        const Eigen::VectorBlock<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1>, Eigen::Dynamic> force_i =
                            force.segment(nc, nc_i);
                        m_i.model.updateForce(d_i, force_i);
                        const pinocchio::JointIndex joint =
                            state_->get_robot()->frames[d_i.frame()].parent;
                        data.fext[joint] = d_i.fext();
                    }
                    else
                    {
                        m_i.model.setZeroForce(d_i);
                    }
                    nc += nc_i;
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
                        const int nc_i = m_i.model.nc();
                        const Eigen::VectorBlock<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, 1>, Eigen::Dynamic> force_i =
                            force.segment(nc, nc_i);
                        m_i.model.updateForce(d_i, force_i);
                        const pinocchio::JointIndex joint =
                            state_->get_robot()->frames[d_i.frame()].parent;
                        data.fext[joint] = d_i.fext();
                        nc += nc_i;
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
            int nc = 0;
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
                    const int nc_i = m_i.model.nc();
                    if (m_i.active)
                    {
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX>> df_dx_i =
                            df_dx.block(nc, 0, nc_i, PS::NDX);
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NU>> df_du_i =
                            df_du.block(nc, 0, nc_i, PS::NU);
                        m_i.model.updateForceDiff(d_i, df_dx_i, df_du_i);
                    }
                    else
                    {
                        m_i.model.setZeroForceDiff(d_i);
                    }
                    nc += nc_i;
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
                        const int nc_i = m_i.model.nc();
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX>> df_dx_i =
                            df_dx.block(nc, 0, nc_i, PS::NDX);
                        const Eigen::Block<const Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::NU>> df_du_i =
                            df_du.block(nc, 0, nc_i, PS::NU);
                        m_i.model.updateForceDiff(d_i, df_dx_i, df_du_i);
                        nc += nc_i;
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
                    switch (m_i.model.type_impl())
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

        const ModelContainer_t &getContacts() const
        {
            return contacts_;
        }

        const State_t &getState() const
        {
            return *state_;
        }

        int nc() const
        {
            return nc_;
        }

        int nc_total() const
        {
            return nc_total_;
        }

    protected:
        State_t *state_;
        ModelContainer_t contacts_;

        int nc_;
        int nc_total_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;
        bool compute_all_contacts_;

    }; // class ContactModelManagerTpl

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_manager_hpp__
