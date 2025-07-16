#ifndef __galileo_core_costs_cost_manager_hpp__
#define __galileo_core_costs_cost_manager_hpp__

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "galileo/core/costs/cost-generic.hpp"
#include "galileo/core/costs/fwd.hpp"

// Despite its name, this file only pertains to cost residuals,
// not just any general costs derived from CostModelBase/CostDataBase

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    struct CostManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    struct CostItemTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using NumScalar = typename PS::NumScalar;

        CostItemTpl(const std::string &name_, const Model_t &model_,
                    const NumScalar &weight_, bool active_ = true)
            : name(name_), model(model_), weight(weight_), active(active_) {}

        std::string name;
        Model_t model;
        NumScalar weight;
        bool active;
    };

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    struct traits<CostManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = CostCollectionTpl<PS>;
        using ModelManager_t = CostModelManagerTpl<PS, CostCollectionTpl>;
        using DataManager_t = CostDataManagerTpl<PS, CostCollectionTpl>;

        using Meta_t = CostTpl<PS, CostCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Item_t = CostItemTpl<PS, CostCollectionTpl>;

        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;

        using L_t = typename PS::VarScalar;
        using Lx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, 1, PS::Options>;
        using Lu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNU_t::Value, 1, PS::Options>;
        using Lxx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Lxu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNDX_t::Value, PS::DimNU_t::Value, PS::Options>;
        using Luu_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNU_t::Value, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    struct traits<CostDataManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    struct traits<CostModelManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    class CostDataManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GALILEO_COST_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        CostDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : L(0.),
              Lx(model_manager.get_ps().ndx_dim.value()),
              Lu(model_manager.get_ps().nu_dim.value()),
              Lxx(model_manager.get_ps().ndx_dim.value(), model_manager.get_ps().ndx_dim.value()),
              Lxu(model_manager.get_ps().ndx_dim.value(), model_manager.get_ps().nu_dim.value()),
              Luu(model_manager.get_ps().nu_dim.value(), model_manager.get_ps().nu_dim.value())
        {
            Lx.setZero();
            Lu.setZero();
            Lxx.setZero();
            Lxu.setZero();
            Luu.setZero();
            for (typename ModelManager_t::ModelContainer_t::const_iterator
                     it = model_manager.get_costs().begin();
                 it != model_manager.get_costs().end(); ++it)
            {
                const Item_t &item = it->second;
                costs.insert(std::make_pair(item.name, item.model.createData(collector)));
            }
        }

        DataContainer_t costs;
        L_t L;
        Lx_t Lx;
        Lu_t Lu;
        Lxx_t Lxx;
        Lxu_t Lxu;
        Luu_t Luu;

    }; // class CostDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class CostCollectionTpl>
    class CostModelManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using NumScalar = typename PS::NumScalar;
        using VarScalar = typename PS::VarScalar;

        CostModelManagerTpl(const PS &ps)
            : ps_(ps), nr_active_dim_(0), nr_total_dim_(0)
        {
        }

        void addCost(const std::string &name, const Model_t &model, const NumScalar &weight,
                     const bool active = true)
        {
            std::pair<typename ModelContainer_t::iterator, bool> ret =
                costs_.insert(std::make_pair(
                    name, Item_t(name, model, weight, active)));
            if (ret.second == false)
            {
                std::cerr << "Warning: we couldn't add the " << name
                          << " cost item, it already existed." << std::endl;
            }
            else if (active)
            {
                nr_active_dim_ += model.get_nr_dim();
                nr_total_dim_ += model.get_nr_dim();
                active_set_.insert(name);
            }
            else if (!active)
            {
                nr_total_dim_ += model.get_nr_dim();
                inactive_set_.insert(name);
            }
        }

        void removeCost(const std::string &name)
        {
            typename ModelContainer_t::iterator it = costs_.find(name);
            if (it != costs_.end())
            {
                nr_active_dim_ -= it->second.model.get_nr_dim();
                nr_total_dim_ -= it->second.model.get_nr_dim();
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
            typename ModelContainer_t::iterator it = costs_.find(name);
            if (it != costs_.end())
            {
                if (active && !it->second.active)
                {
                    nr_active_dim_ += it->second.model.get_nr_dim();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    nr_active_dim_ -= it->second.model.get_nr_dim();
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
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.L = VarScalar(0.);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                end_d = data.costs.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived(), u.derived());
                    data.L += m_i.weight * d_i.L();
                }
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            data.L = VarScalar(0.);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                end_d = data.costs.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived());
                    data.L += m_i.weight * d_i.L();
                }
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.Lx.setZero();
            data.Lu.setZero();
            data.Lxx.setZero();
            data.Lxu.setZero();
            data.Luu.setZero();

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                end_d = data.costs.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived(), u.derived());
                    data.Lx += m_i.weight * d_i.Lx();
                    data.Lu += m_i.weight * d_i.Lu();
                    data.Lxx += m_i.weight * d_i.Lxx();
                    data.Lxu += m_i.weight * d_i.Lxu();
                    data.Luu += m_i.weight * d_i.Luu();
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            data.Lx.setZero();
            data.Lxx.setZero();

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = costs_.begin(), end_m = costs_.end(), it_d = data.costs.begin(),
                end_d = data.costs.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived());
                    data.Lx += m_i.weight * d_i.Lx();
                    data.Lxx += m_i.weight * d_i.Lxx();
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

        const ModelContainer_t &get_costs() const
        {
            return costs_;
        }

        const std::set<std::string> &get_active_set() const
        {
            return active_set_;
        }

        const std::set<std::string> &get_inactive_set() const
        {
            return inactive_set_;
        }

        const bool get_cost_status(const std::string &name) const
        {
            typename ModelContainer_t::const_iterator it =
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

        const int get_nr_active() const
        {
            return nr_active_dim_.value();
        }

        const int get_nr_total() const
        {
            return nr_total_dim_.value();
        }

    protected:
        const PS &ps_;
        ModelContainer_t costs_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

        DimensionTpl<Eigen::Dynamic> nr_active_dim_;
        DimensionTpl<Eigen::Dynamic> nr_total_dim_;

    }; // class CostModelManagerTpl

} // namespace galileo

#endif // __galileo_core_costs_cost_manager_hpp__
