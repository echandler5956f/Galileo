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

#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace multibody
    {

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        struct ContactItemTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactCollection = ContactCollectionTpl<PS>;

            using ContactModel = ContactModelTpl<PS, ContactCollectionTpl>;
            using ContactData = ContactDataTpl<PS, ContactCollectionTpl>;

            ContactItemTpl() {}
            ContactItemTpl(const std::string &name, const ContactModel &contact, bool active = true)
                : name(name), contact(contact), active(active) {}

            std::string name;
            ContactModel contact;
            bool active;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        class ContactDataManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactCollection = ContactCollectionTpl<PS>;

            using ContactData = ContactDataTpl<PS, ContactCollectionTpl>;
            using ContactDataContainer = std::map<std::string, ContactData>;

            using ContactDerived = typename traits<ContactData>::ContactDerived;

            GALILEO_FORCE_DATA_TYPEDEF(ContactDerived);
            GALILEO_CONTACT_DATA_TYPEDEF(ContactDerived);

            Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NV, PS::Options> Jc;
            Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, 1, PS::Options> a0;
            Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX, PS::Options> da0_dx;
            Eigen::Matrix<typename PS::VarScalar, PS::NV, 1, PS::Options> dv;
            Eigen::Matrix<typename PS::VarScalar, PS::NV, PS::NDX, PS::Options> ddv_dx;

            ContactDataContainer contacts;
            using fext = GALILEO_ALIGNED_STD_VECTOR(Force_t);

        }; // class ContactDataManagerTpl

        template <typename PhaseSpec,
                  template <typename PS> class ContactCollectionTpl>
        class ContactModelManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactCollection = ContactCollectionTpl<PS>;

            using ContactData = ContactDataTpl<PS, ContactCollectionTpl>;
            using ContactModel = ContactModelTpl<PS, ContactCollectionTpl>;

            using ContactItem = ContactItemTpl<PS, ContactCollectionTpl>;

            using ContactModelContainer = std::map<std::string, ContactItem>;
            using ContactDataContainer = std::map<std::string, ContactData>;

            using ContactDataManager = ContactDataManagerTpl<PS, ContactCollectionTpl>;

            using RobotData = typename PS::RobotData_t;

            using ForceIterator = typename galileo::container::aligned_vector<typename PS::Force_t>::iterator

            ContactModelManagerTpl() {}

            void addContact(const std::string &name, const ContactModel &contact, bool active = true)
            {
                std::pair<typename ContactModelContainer::iterator, bool> ret =
                    contacts_.insert(std::make_pair(
                        name, ContactItem(name, contact, active)));
                if (ret.second == false)
                {
                    std::cout << "Warning: we couldn't add the " << name
                              << " contact item, it already existed." << std::endl;
                }
                else if (active)
                {
                    nc_ += contact.nc();
                    nc_total_ += contact.nc();
                    active_set_.insert(name);
                }
                else if (!active)
                {
                    nc_total_ += contact.nc();
                    inactive_set_.insert(name);
                }
            }

            void removeContact(const std::string &name)
            {
                typename ContactModelContainer::iterator it = contacts_.find(name);
                if (it != contacts_.end())
                {
                    nc_ -= it->second.contact.nc();
                    nc_total_ -= it->second.contact.nc();
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
                typename ContactModelContainer::iterator it = contacts_.find(name);
                if (it != contacts_.end())
                {
                    if (active && !it->second.active)
                    {
                        nc_ += it->second.contact.nc();
                        active_set_.insert(name);
                        inactive_set_.erase(name);
                        it->second.active = active;
                    }
                    else if (!active && it->second.active)
                    {
                        nc_ -= it->second.contact.nc();
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
            void calc(ContactDataManager &data, const Eigen::MatrixBase<StateVectorType> &x)
            {
                int nc = 0;
                typename ContactModelContainer::iterator it_m, end_m;
                typename ContactDataContainer::iterator it_d, end_d;
                if (compute_all_contacts_)
                {
                    for (it_m = contacts_.begin(), end_m = contacts_.end(),
                        it_d = data.contacts.begin(), end_d = data.contacts.end();
                         it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                    {
                        ContactItem &m_i = it_m->second;
                        const int nc_i = m_i.contact.nc();
                        if (m_i.active)
                        {
                            ContactData &d_i = it_d->second;

                            m_i.contact.calc(d_i, x.derived());
                            data.a0().segment(nc, nc_i) = d_i.a0();
                            data.Jc().block(nc, 0, nc_i, PS::NV) = d_i.Jc();
                        }
                        else
                        {
                            data.a0().segment(nc, nc_i).setZero();
                            data.Jc().block(nc, 0, nc_i, PS::NV).setZero();
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
                        ContactItem &m_i = it_m->second;
                        if (m_i.active)
                        {
                            ContactData &d_i = it_d->second;

                            m_i.contact.calc(d_i, x.derived());
                            const int nc_i = m_i.contact.nc();
                            data.a0().segment(nc, nc_i) = d_i.a0();
                            data.Jc().block(nc, 0, nc_i, PS::NV) = d_i.Jc();
                            nc += nc_i;
                        }
                    }
                }
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataManager &data, const Eigen::MatrixBase<StateVectorType> &x)
            {
                int nc = 0;
                typename ContactModelContainer::iterator it_m, end_m;
                typename ContactDataContainer::iterator it_d, end_d;
                if (compute_all_contacts_)
                {
                    for (it_m = contacts_.begin(), end_m = contacts_.end(),
                        it_d = data.contacts.begin(), end_d = data.contacts.end();
                         it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                    {
                        ContactItem &m_i = it_m->second;
                        const int nc_i = m_i.contact.nc();
                        if (m_i.active)
                        {
                            ContactData &d_i = it_d->second;

                            m_i.contact.calcDiff(d_i, x.derived());
                            data.da0_dx().block(nc, 0, nc_i, PS::NDX) = d_i.da0_dx();
                        }
                        else
                        {
                            data.da0_dx().block(nc, 0, nc_i, PS::NDX).setZero();
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
                        ContactItem &m_i = it_m->second;
                        if (m_i.active)
                        {
                            ContactData &d_i = it_d->second;

                            m_i.contact.calcDiff(d_i, x.derived());
                            const int nc_i = m_i.contact.nc();
                            data.da0_dx().block(nc, 0, nc_i, PS::NDX) = d_i.da0_dx();
                            nc += nc_i;
                        }
                    }
                }
            }

            template <typename VectorNvType>
            void updateAcceleration(ContactDataManager &data, const Eigen::MatrixBase<VectorNvType> &dv) const
            {
                data.dv() = dv.derived();
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataManager &data, const Eigen::MatrixBase<ForceVectorType> &force)
            {
                for (ForceIterator it = data.fext.begin(); it != data.fext.end(); ++it)
                {
                    *it = typename PS::Force_t::Zero();
                }

                std::size_t nc = 0;
                typename ContactModelContainer::iterator it_m, end_m;
                typename ContactDataContainer::iterator it_d, end_d;
                if (compute_all_contacts_)
                {
                    for (it_m = contacts_.begin(), end_m = contacts_.end(),
                        it_d = data.contacts.begin(), end_d = data.contacts.end();
                         it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                    {
                        ContactItem &m_i = it_m->second;
                        ContactData &d_i = it_d->second;
                        const int nc_i = m_i.contact.nc();
                        if (m_i.active)
                        {
                            const Eigen::VectorBlock<const VectorXs, Eigen::Dynamic> force_i =
                                force.segment(nc, nc_i);
                            m_i.contact.updateForce(d_i, force_i);
                            const pinocchio::JointIndex joint =
                                state_->get_pinocchio()->frames[d_i.frame].parent;
                            data.fext[joint] = d_i.fext;
                        }
                        else
                        {
                            m_i.contact.setZeroForce(d_i);
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
                        ContactItem &m_i = it_m->second;
                        ContactData &d_i = it_d->second;
                        if (m_i.active)
                        {
                            const int nc_i = m_i.contact.nc();
                            const Eigen::VectorBlock<const VectorXs, Eigen::Dynamic> force_i =
                                force.segment(nc, nc_i);
                            m_i.contact.updateForce(d_i, force_i);
                            const pinocchio::JointIndex joint =
                                state_->get_pinocchio()->frames[d_i.frame].parent;
                            data.fext[joint] = d_i.fext;
                            nc += nc_i;
                        }
                        else
                        {
                            m_i.contact.setZeroForce(d_i);
                        }
                    }
                }
            }

            template <typename MatrixNvNdxType>
            void updateAccelerationDiff(ContactDataManager &data, const Eigen::MatrixBase<MatrixNvNdxType> &ddv_dx) const
            {
                data.ddv_dx() = ddv_dx.derived();
            }

            template <typename MatrixNcNdxType, typename MatrixNcNduType>
            void updateForceDiff(ContactDataManager &data, const Eigen::MatrixBase<MatrixNcNdxType> &df_dx, const Eigen::MatrixBase<MatrixNcNduType> &df_du) const
            {
                int nc = 0;
                typename ContactModelContainer::const_iterator it_m, end_m;
                typename ContactDataContainer::const_iterator it_d, end_d;
                if (compute_all_contacts_)
                {
                    for (it_m = contacts_.begin(), end_m = contacts_.end(),
                        it_d = data.contacts.begin(), end_d = data.contacts.end();
                         it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                    {
                        ContactItem &m_i = it_m->second;
                        ContactData &d_i = it_d->second;
                        const int nc_i = m_i.contact.nc();
                        if (m_i.active)
                        {
                            const Eigen::Block<const Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX>> df_dx_i =
                                df_dx.block(nc, 0, nc_i, PS::NDX);
                            const Eigen::Block<const Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NU>> df_du_i =
                                df_du.block(nc, 0, nc_i, PS::NU);
                            m_i.contact.updateForceDiff(d_i, df_dx_i, df_du_i);
                        }
                        else
                        {
                            m_i.contact.setZeroForceDiff(d_i);
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
                        ContactItem &m_i = it_m->second;
                        ContactData &d_i = it_d->second;
                        if (m_i.active)
                        {
                            const int nc_i = m_i.contact.nc();
                            const Eigen::Block<const Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX>> df_dx_i =
                                df_dx.block(nc, 0, nc_i, PS::NDX);
                            const Eigen::Block<const Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NU>> df_du_i =
                                df_du.block(nc, 0, nc_i, PS::NU);
                            m_i.contact.updateForceDiff(d_i, df_dx_i, df_du_i);
                            nc += nc_i;
                        }
                        else
                        {
                            m_i.contact.setZeroForceDiff(d_i);
                        }
                    }
                }
            }

            void updateRneaDiff(ContactDataManager &data, RobotData &robot_data) const
            {
                typename ContactModelContainer::const_iterator it_m, end_m;
                typename ContactDataContainer::const_iterator it_d, end_d;
                for (it_m = contacts_.begin(), end_m = contacts_.end(),
                    it_d = data.contacts.begin(), end_d = data.contacts.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    ContactItem &m_i = it_m->second;
                    ContactData &d_i = it_d->second;
                    if (m_i.active)
                    {
                        switch (m_i.contact.type())
                        {
                        case pinocchio::ReferenceFrame::LOCAL:
                            break;
                        case pinocchio::ReferenceFrame::WORLD:
                        case pinocchio::ReferenceFrame::LOCAL_WORLD_ALIGNED:
                            data.dtau_dq() += d_i.dtau_dq();
                            break;
                        }
                    }
                }
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
            using State = typename PS::State_t;
            State *state_;
            ContactModelContainer contacts_;

            int nc_;
            int nc_total_;

            std::set<std::string> active_set_;
            std::set<std::string> inactive_set_;
            bool compute_all_contacts_;

        }; // class ContactModelManagerTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_manager_hpp__