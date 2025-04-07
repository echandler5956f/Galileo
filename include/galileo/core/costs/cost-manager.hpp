#ifndef __galileo_core_costs_cost_manager_hpp__
#define __galileo_core_costs_cost_manager_hpp__

#include <iostream>
#include <string>
#include <map>

#include "galileo/core/costs/fwd.hpp"
#include "galileo/core/costs/cost-generic.hpp"

// Despite its name, this file only pertains to cost residuals,
// not just any general costs derived from CostModelBase/CostDataBase

namespace galileo
{
    namespace core
    {

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        struct CostItemTpl
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostCollection = CostCollectionTpl<PS>;

            using CostModel = CostModelTpl<PS, CostCollectionTpl>;
            using CostData = CostDataTpl<PS, CostCollectionTpl>;

            CostItemTpl() {}
            CostItemTpl(const std::string &name, const CostModel &cost, const PS::NumScalar &weight, bool active = true)
                : name(name), cost(cost), weight(weight), active(active) {}

            std::string name;
            CostModel cost;
            PS::NumScalar weight;
            bool active;
        };

        template <typename PhaseSpec,
                  template <typename PS> class CostCollectionTpl>
        class CostDataManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostCollection = CostCollectionTpl<PS>;

            using CostData = CostDataTpl<PS, CostCollectionTpl>;
            using CostDataContainer = std::map<std::string, CostData>;

            using CostDerived = typename traits<CostData>::CostDerived;

            GALILEO_COST_DATA_TYPEDEF(CostDerived);

            CostDataContainer costs;
            Eigen::Map<L_t> L;
            Eigen::Map<Lx_t> Lx;
            Eigen::Map<Lu_t> Lu;
            Eigen::Map<Lxx_t> Lxx;
            Eigen::Map<Lxu_t> Lxu;
            Eigen::Map<Luu_t> Luu;

        }; // class CostDataManagerTpl

        template <typename PS,
                  template <typename PS> class CostCollectionTpl>
        class CostModelManagerTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using CostCollection = CostCollectionTpl<PS>;

            using CostGeneric = CostTpl<PS, CostCollectionTpl>;
            using CostModel = CostModelTpl<PS, CostCollectionTpl>;
            using CostData = CostDataTpl<PS, CostCollectionTpl>;

            using CostItem = CostItemTpl<PS, CostCollectionTpl>;

            using CostModelContainer = std::map<std::string, CostItem>;
            using CostDataContainer = std::map<std::string, CostData>;

            using CostModelManager = CostModelManagerTpl<PS, CostCollectionTpl>;
            using CostDataManager = CostDataManagerTpl<PS, CostCollectionTpl>;

            CostModelManagerTpl() {}

            void addCost(const std::string &name, const CostModel &cost, const PS::NumScalar &weight, const bool active = true)
            {
                std::pair<typename CostModelContainer::iterator, bool> ret =
                    costs_.insert(std::make_pair(
                        name, CostItem(name, cost, weight, active)));
                if (ret.second == false)
                {
                    std::cerr << "Warning: we couldn't add the " << name
                              << " cost item, it already existed." << std::endl;
                }
                else if (active)
                {
                    nr_ += cost.nr();
                    nr_total_ += cost.nr();
                    active_set_.insert(name);
                }
                else if (!active)
                {
                    nr_total_ += cost.nr();
                    inactive_set_.insert(name);
                }
            }

            void removeCost(const std::string &name)
            {
                typename CostModelContainer::iterator it = costs_.find(name);
                if (it != costs_.end())
                {
                    nr_ -= it->second.cost.nr();
                    nr_total_ -= it->second.cost.nr();
                    costs_.erase(it);
                    active_set_.erase(name);
                    inactive_set_.erase(name);
                }
                else
                {
                    std::cerr << "Warning: we couldn't remove the " << name
                              << " cost item, it doesn't exist." << std::endl;
                }
            }

            void changeCostStatus(const std::string &name, bool active)
            {
                typename CostModelContainer::iterator it = costs_.find(name);
                if (it != costs_.end())
                {
                    if (active && !it->second.active)
                    {
                        nr_ += it->second.cost.nr();
                        active_set_.insert(name);
                        inactive_set_.erase(name);
                        it->second.active = active;
                    }
                    else if (!active && it->second.active)
                    {
                        nr_ -= it->second.cost.nr();
                        active_set_.erase(name);
                        inactive_set_.insert(name);
                        it->second.active = active;
                    }
                }
                else
                {
                    std::cerr << "Warning: we couldn't change the status of the " << name
                              << " cost item, it doesn't exist." << std::endl;
                }
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calc(CostDataManager &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u)
            {
                data.cost = typename PS::NumScalar(0.);

                typename CostModelContainer::iterator it_m, end_m;
                typename CostDataContainer::iterator it_d, end_d;
                for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                    end_d = data.costs.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const CostItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        CostData &d_i = it_d->second;

                        m_i.cost.calc(d_i, x.derived(), u.derived());
                        data.L() += m_i.weight * d_i.L();
                    }
                }
            }

            template <typename StateVectorType>
            void calc(CostDataManager &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
            {
                data.L() = typename PS::NumScalar(0.);

                typename CostModelContainer::iterator it_m, end_m;
                typename CostDataContainer::iterator it_d, end_d;
                for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                    end_d = data.costs.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const CostItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        CostData &d_i = it_d->second;

                        m_i.cost.calc(d_i, x.derived());
                        data.L() += m_i.weight * d_i.L();
                    }
                }
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(CostDataManager &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u)
            {
                data.Lx().setZero();
                data.Lu().setZero();
                data.Lxx().setZero();
                data.Lxu().setZero();
                data.Luu().setZero();

                typename CostModelContainer::iterator it_m, end_m;
                typename CostDataContainer::iterator it_d, end_d;
                for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                    end_d = data.costs.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const CostItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        CostData &d_i = it_d->second;

                        m_i.cost.calcDiff(d_i, x.derived(), u.derived());
                        data.Lx() += m_i.weight * d_i.Lx();
                        data.Lu() += m_i.weight * d_i.Lu();
                        data.Lxx() += m_i.weight * d_i.Lxx();
                        data.Lxu() += m_i.weight * d_i.Lxu();
                        data.Luu() += m_i.weight * d_i.Luu();
                    }
                }
            }

            template <typename StateVectorType>
            void calcDiff(CostDataManager &data,
                          const Eigen::MatrixBase<StateVectorType> &x)
            {
                data.Lx().setZero();
                data.Lxx().setZero();

                typename CostModelContainer::iterator it_m, end_m;
                typename CostDataContainer::iterator it_d, end_d;
                for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                    end_d = data.costs.end();
                     it_m != end_m || it_d != end_d; ++it_m, ++it_d)
                {
                    const CostItem &m_i = it_m->second;
                    if (m_i.active)
                    {
                        CostData &d_i = it_d->second;

                        m_i.cost.calcDiff(d_i, x.derived());
                        data.Lx() += m_i.weight * d_i.Lx();
                        data.Lxx() += m_i.weight * d_i.Lxx();
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

            bool getCostStatus(const std::string &name) const
            {
                typename CostModelContainer::const_iterator it =
                    costs_.find(name);
                if (it != costs_.end())
                {
                    return it->second.active;
                }
                else
                {
                    std::cout << "Warning: we couldn't get the status of the " << name
                              << " cost item, it doesn't exist." << std::endl;
                    return false;
                }
            }

        protected:
            CostModelContainer costs_;

            std::size_t nr_;
            std::size_t nr_total_;

            std::set<std::string> active_set_;
            std::set<std::string> inactive_set_;

        }; // class CostModelManagerTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_manager_hpp__
