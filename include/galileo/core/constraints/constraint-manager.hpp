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

        template <typename PhaseSpec,
                  template <typename PS> class ConstraintCollectionTpl>
        struct ConstraintItemTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintCollection = ConstraintCollectionTpl<PS>;

            using ConstraintModel = ConstraintModelTpl<PS, ConstraintCollectionTpl>;
            using ConstraintData = ConstraintDataTpl<PS, ConstraintCollectionTpl>;

            ConstraintItemTpl() {}
            ConstraintItemTpl(const std::string &name, const ConstraintModel &constraint, bool active = true)
                : name(name), constraint(constraint), active(active) {}

            std::string name;
            ConstraintModel constraint;
            bool active;
        };

        template <typename PhaseSpec,
                  template <typename PS> class ConstraintCollectionTpl>
        class ConstraintDataManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintCollection = ConstraintCollectionTpl<PS>;

            using ConstraintData = ConstraintDataTpl<PS, ConstraintCollectionTpl>;
            using ConstraintDataContainer = std::map<std::string, ConstraintData>;

            using ConstraintDerived = typename traits<ConstraintData>::ConstraintDerived;

            GALILEO_CONSTRAINT_DATA_TYPEDEF(ConstraintDerived);

            ConstraintDataContainer constraints;
            Eigen::Map<H_t> H;
            Eigen::Map<Hx_t> Hx;
            Eigen::Map<Hu_t> Hu;
            Eigen::Map<G_t> G;
            Eigen::Map<Gx_t> Gx;
            Eigen::Map<Gu_t> Gu;

        }; // class ConstraintDataManagerTpl

        template <typename PS,
                  template <typename PS> class ConstraintCollectionTpl>
        class ConstraintModelManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintCollection = ConstraintCollectionTpl<PS>;

            using ConstraintGeneric = ConstraintTpl<PS, ConstraintCollectionTpl>;
            using ConstraintModel = ConstraintModelTpl<PS, ConstraintCollectionTpl>;
            using ConstraintData = ConstraintDataTpl<PS, ConstraintCollectionTpl>;

            using ConstraintItem = ConstraintItemTpl<PS, ConstraintCollectionTpl>;

            using ConstraintModelContainer = std::map<std::string, ConstraintItem>;
            using ConstraintDataContainer = std::map<std::string, ConstraintData>;

            using ConstraintModelManager = ConstraintModelManagerTpl<PS, ConstraintCollectionTpl>;
            using ConstraintDataManager = ConstraintDataManagerTpl<PS, ConstraintCollectionTpl>;

            using BoundVector_t = typename traits<ConstraintGeneric>::BoundVector_t;

            ConstraintModelManagerTpl() {}

            void add_constraint(const std::string &name, const ConstraintModel &constraint, const bool active = true)
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
                    ng_ += constraint.ng();
                    nh_ += constraint.nh();
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
                    ng_ -= it->second.constraint.ng();
                    nh_ -= it->second.constraint.nh();
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
                        ng_ += it->second.constraint.ng();
                        nh_ += it->second.constraint.nh();
                        active_set_.insert(name);
                        inactive_set_.erase(name);
                        it->second.active = active;
                        lb_.resize(ng_);
                        ub_.resize(ng_);
                    }
                    else if (!active && it->second.active)
                    {
                        ng_ -= it->second.constraint.ng();
                        nh_ -= it->second.constraint.nh();
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
            void calc(ConstraintDataManager &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u)
            {
                std::size_t ng_i = 0;
                std::size_t nh_i = 0;

                typename ConstraintModelContainer::iterator it_m, end_m;
                typename ConstraintDataContainer::iterator it_d, end_d;
                for (it_m = constraints_.begin(), end_m = constraints_.end(),
                    it_d = data.constraints.begin(), end_d = data.constraints.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const ConstraintItem &m_i = it_m->second;
                    if (m_i->active)
                    {
                        ConstraintData &d_i = it_d->second;

                        m_i.constraint.calc(d_i, x.derived(), u.derived());
                        const std::size_t ng = m_i.constraint.ng();
                        const std::size_t nh = m_i.constraint.nh();
                        data.G.segment(ng_i, ng) = d_i.G();
                        data.H.segment(nh_i, nh) = d_i.H();
                        lb_.segment(ng_i, ng) = m_i.constraint.lb();
                        ub_.segment(ng_i, ng) = m_i.constraint.ub();
                        ng_i += ng;
                        nh_i += nh;
                    }
                }
            }

            template <typename StateVectorType>
            void calc(ConstraintDataManager &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                std::size_t ng_i = 0;
                std::size_t nh_i = 0;

                typename ConstraintModelContainer::iterator it_m, end_m;
                typename ConstraintDataContainer::iterator it_d, end_d;
                for (it_m = constraints_.begin(), end_m = constraints_.end(),
                    it_d = data.constraints.begin(), end_d = data.constraints.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const ConstraintItem &m_i = it_m->second;
                    if (m_i->active)
                    {
                        ConstraintData &d_i = it_d->second;

                        m_i.constraint.calc(d_i, x.derived());
                        const std::size_t ng = m_i.constraint.ng();
                        const std::size_t nh = m_i.constraint.nh();
                        data.G.segment(ng_i, ng) = d_i.G();
                        data.H.segment(nh_i, nh) = d_i.H();
                        lb_.segment(ng_i, ng) = m_i.constraint.lb();
                        ub_.segment(ng_i, ng) = m_i.constraint.ub();
                        ng_i += ng;
                        nh_i += nh;
                    }
                }
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataManager &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u)
            {
                std::size_t ng_i = 0;
                std::size_t nh_i = 0;

                typename ConstraintModelContainer::iterator it_m, end_m;
                typename ConstraintDataContainer::iterator it_d, end_d;
                for (it_m = constraints_.begin(), end_m = constraints_.end(),
                    it_d = data.constraints.begin(), end_d = data.constraints.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    ConstraintItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        ConstraintData &d_i = it_d->second;

                        m_i.constraint.calcDiff(d_i, x.derived(), u.derived());
                        const std::size_t ng = m_i.constraint.ng();
                        const std::size_t nh = m_i.constraint.nh();
                        data.Gx.block(ng_i, 0, ng, PS::NDX) = d_i.Gx();
                        data.Gu.block(ng_i, 0, ng, PS::NU) = d_i.Gu();
                        data.Hx.block(nh_i, 0, nh, PS::NDX) = d_i.Hx();
                        data.Hu.block(nh_i, 0, nh, PS::NU) = d_i.Hu();
                        ng_i += ng;
                        nh_i += nh;
                    }
                }
            }

            template <typename StateVectorType>
            void calcDiff(ConstraintDataManager &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                std::size_t ng_i = 0;
                std::size_t nh_i = 0;

                typename ConstraintModelContainer::iterator it_m, end_m;
                typename ConstraintDataContainer::iterator it_d, end_d;
                for (it_m = constraints_.begin(), end_m = constraints_.end(),
                    it_d = data.constraints.begin(), end_d = data.constraints.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    ConstraintItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        ConstraintData &d_i = it_d->second;

                        m_i.constraint.calcDiff(d_i, x.derived());
                        const std::size_t ng = m_i.constraint.ng();
                        const std::size_t nh = m_i.constraint.nh();
                        data.Gx.block(ng_i, 0, ng, PS::NDX) = d_i.Gx();
                        data.Hx.block(nh_i, 0, nh, PS::NDX) = d_i.Hx();
                        ng_i += ng;
                        nh_i += nh;
                    }
                }
            }

            const std::set<std::string> &getActiveSet() const
            {
                return active_set_;
            }

            const std::set<std::string> &getInactiveSet() const
            {
                return inactive_set_;
            }

            const BoundVector_t &get_lb() const
            {
                return lb_;
            }

            const BoundVector_t &get_ub() const
            {
                return ub_;
            }

            bool getConstraintStatus(const std::string &name) const
            {
                typename ConstraintModelContainer::const_iterator it =
                    constraints_.find(name);
                if (it != constraints_.end())
                {
                    return it->second.active;
                }
                else
                {
                    std::cout << "Warning: we couldn't get the status of the " << name
                              << " constraint item, it doesn't exist." << std::endl;
                    return false;
                }
            }

        protected:
            ConstraintModelContainer constraints_;
            BoundVector_t lb_;
            BoundVector_t ub_;

            std::size_t nh_;
            std::size_t ng_;

            std::set<std::string> active_set_;
            std::set<std::string> inactive_set_;

        }; // class ConstraintModelManagerTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_manager_hpp__
