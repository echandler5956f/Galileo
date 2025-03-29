#ifndef __galileo_core_constraints_constraint_manager_hpp__
#define __galileo_core_constraints_constraint_manager_hpp__

#include <iostream>
#include <string>
#include <map>

#include "galileo/core/constraints/fwd.hpp"
#include "galileo/core/constraints/constraint-generic.hpp"

namespace galileo
{
    namespace core
    {

        template <typename _VarScalar, typename _NumScalar, int _Options,
                  template <typename V, typename N, int O> class _ConstraintCollectionTpl>
        struct ConstraintItemTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            using Options = _Options;

            using ConstraintCollection = _ConstraintCollectionTpl<VarScalar, NumScalar, Options>;

            using ConstraintModel = ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollection>;
            using ConstraintData = ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollection>;

            ConstraintItemTpl() {}
            ConstraintItemTpl(const std::string &name, const ConstraintModel &constraint, bool active = true)
                : name(name), constraint(constraint), active(active) {}

            std::string name;
            ConstraintModel constraint;
            bool active;
        };

        template <typename _VarScalar, typename _NumScalar, int _Options,
                  template <typename V, typename N, int O> class _ConstraintCollectionTpl>
        class ConstraintModelManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            using Options = _Options;

            using ConstraintCollection = _ConstraintCollectionTpl<VarScalar, NumScalar, Options>;

            using ConstraintModel = ConstraintModelTpl<VarScalar, NumScalar, Options, ConstraintCollection>;
            using ConstraintData = ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollection>;

            using ConstraintItem = ConstraintItemTpl<VarScalar, NumScalar, Options, ConstraintCollection>;

            using ConstraintModelContainer = std::map<std::string, ConstraintItem>;
            using ConstraintDataContainer = std::map<std::string, ConstraintData>;

            ConstraintModelManagerTpl() {}

            void add_constraint(const std::string &name, const ConstraintModel &constraint, bool active = true)
            {
                std::pair<typename ConstraintModelContainer::iterator, bool> ret =
                    constraints_.insert(std::make_pair(
                        name, ConstraintItem(name, constraint, active)));
                if (ret.second == false)
                {
                    std::cout << "Warning: we couldn't add the " << name
                              << " constraint item, it already existed." << std::endl;
                }
                else if (active)
                {
                    ng_ += constraint.get_ng();
                    nh_ += constraint.get_nh();
                    active_set_.insert(name);
                    lb_.resize(ng_);
                    ub_.resize(ng_);
                }
                else if (!active)
                {
                    inactive_set_.insert(name);
                }
            }

            void remove_constraint(const std::string &name)
            {
                typename ConstraintModelContainer::iterator it = constraints_.find(name);
                if (it != constraints_.end())
                {
                    ng_ -= it->second.constraint.get_ng();
                    nh_ -= it->second.constraint.get_nh();
                    constraints_.erase(it);
                    inactive_set_.erase(name);
                    lb_.resize(ng_);
                    ub_.resize(ng_);
                }
                else
                {
                    std::cout << "Warning: we couldn't remove the " << name
                              << " constraint item, it doesn't exist." << std::endl;
                }
            }

            void changeConstraintStatus(const std::string &name, bool active)
            {
                typename ConstraintModelContainer::iterator it = constraints_.find(name);
                if (it != constraints_.end())
                {
                    if (active && !it->second.active)
                    {
                        ng_ += it->second.constraint.get_ng();
                        nh_ += it->second.constraint.get_nh();
                        active_set_.insert(name);
                        inactive_set_.erase(name);
                        it->second.active = active;
                        lb_.resize(ng_);
                        ub_.resize(ng_);
                    }
                    else if (!active && it->second.active)
                    {
                        ng_ -= it->second.constraint.get_ng();
                        nh_ -= it->second.constraint.get_nh();
                        active_set_.erase(name);
                        inactive_set_.insert(name);
                        it->second.active = active;
                        lb_.resize(ng_);
                        ub_.resize(ng_);
                    }
                }
                else
                {
                    std::cout << "Warning: we couldn't change the status of the " << name
                              << " constraint item, it doesn't exist." << std::endl;
                }
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<ControlVectorType> &u)
            {
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(const Eigen::MatrixBase<StateVectorType> &x, const Eigen::MatrixBase<ControlVectorType> &u)
            {
            }

            const std::set<std::string> &getActiveSet() const
            {
                return active_set_;
            }

            const std::set<std::string> &getInactiveSet() const
            {
                return inactive_set_;
            }

            template <typename VectorType>
            Eigen::MatrixBase<VectorType> get_lb() const
            {
                return lb_;
            }

            template <typename VectorType>
            Eigen::MatrixBase<VectorType> get_ub() const
            {
                return ub_;
            }

            bool getConstraintStatus(const std::string &name) const
            {
                return constraints_.at(name).active;
            }

        protected:
            ConstraintModelContainer constraints_;
            Eigen::Matrix<NumScalar, Eigen::Dynamic, 1> lb_;
            Eigen::Matrix<NumScalar, Eigen::Dynamic, 1> ub_;

            std::size_t nh_;
            std::size_t ng_;

            std::set<std::string> active_set_;
            std::set<std::string> inactive_set_;

        }; // class ConstraintModelManagerTpl

        template <typename _VarScalar, typename _NumScalar, int _Options,
                  template <typename V, typename N, int O> class _ConstraintCollectionTpl>
        class ConstraintDataManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            using Options = _Options;

            using ConstraintCollection = _ConstraintCollectionTpl<VarScalar, NumScalar, Options>;

            using ConstraintData = ConstraintDataTpl<VarScalar, NumScalar, Options, ConstraintCollection>;
            using ConstraintDataContainer = std::map<std::string, ConstraintData>;

        }; // class ConstraintDataManagerTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_manager_hpp__
