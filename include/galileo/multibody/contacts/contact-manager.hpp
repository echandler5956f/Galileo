#ifndef __galileo_multibody_contacts_contact_manager_hpp__
#define __galileo_multibody_contacts_contact_manager_hpp__

#include <iostream>
#include <string>
#include <map>

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/contacts/contact-generic.hpp"

namespace galileo
{

    namespace multibody
    {

        template <typename PhaseSpec>
        struct ContactItemTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactModel = ContactModelTpl<PS>;
            using ContactData = ContactDataTpl<PS>;

            ContactItemTpl() {}
            ContactItemTpl(const std::string &name, const ContactModel &contact, bool active = true)
                : name(name), contact(contact), active(active) {}

            std::string name;
            ContactModel contact;
            bool active;
        };

        template <typename PhaseSpec>
        class ContactModelManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactModel = ContactModelTpl<PS>;
            using ContactData = ContactDataTpl<PS>;

            using ContactItem = ContactItemTpl<PS>;

            using ContactModelContainer = std::map<std::string, ContactItem>;
            using ContactDataContainer = std::map<std::string, ContactData>;

            using ContactDataManager = ContactDataManagerTpl<PS>;

            using RobotData = typename PS::RobotData_t;

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
                    nc_ += contact.get_nc();
                    active_set_.insert(name);
                }
                else if (!active)
                {
                    inactive_set_.insert(name);
                }
            }

            void removeContact(const std::string &name)
            {
                typename ContactModelContainer::iterator it = contacts_.find(name);
                if (it != contacts_.end())
                {
                    nc_ -= it->second.contact.get_nc();
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
                        nc_ += it->second.contact.get_nc();
                        active_set_.insert(name);
                        inactive_set_.erase(name);
                        it->second.active = active;
                    }
                    else if (!active && it->second.active)
                    {
                        nc_ -= it->second.contact.get_nc();
                        active_set_.erase(name);
                        inactive_set_.insert(name);
                        it->second.active = active;
                    }
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
            }

            template <typename StateVectorType>
            void calcDiff(ContactDataManager &data, const Eigen::MatrixBase<StateVectorType> &x)
            {
            }

            template <typename VectorNvType>
            void updateAcceleration(ContactDataManager &data, const Eigen::MatrixBase<VectorNvType> &dv) const
            {
            }

            template <typename ForceVectorType>
            void updateForce(ContactDataManager &data, const Eigen::MatrixBase<ForceVectorType> &f)
            {
            }

            template <typename MatrixNvNdxType>
            void updateAccelerationDiff(ContactDataManager &data, const Eigen::MatrixBase<MatrixNvNdxType> &ddv_dx) const
            {
            }

            template <typename MatrixNcNdxType, typename MatrixNcNduType>
            void updateForceDiff(ContactDataManager &data, const Eigen::MatrixBase<MatrixNcNdxType> &df_dx, const Eigen::MatrixBase<MatrixNcNduType> &df_du) const
            {
            }

            void updateRneaDiff(ContactDataManager &data, RobotData &robot_data) const
            {
            }

        protected:
            ContactModelContainer contacts_;

            std::size_t nc_;

            std::set<std::string> active_set_;
            std::set<std::string> inactive_set_;

        }; // class ContactModelManagerTpl

        template <typename PhaseSpec>
        class ContactDataManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ContactModel = ContactModelTpl<PS>;
            using ContactData = ContactDataTpl<PS>;

            using ContactItem = ContactItemTpl<PS>;

            using ContactModelContainer = std::map<std::string, ContactItem>;
            using ContactDataContainer = std::map<std::string, ContactData>;

        }; // class ContactDataManagerTpl

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_contacts_contact_manager_hpp__