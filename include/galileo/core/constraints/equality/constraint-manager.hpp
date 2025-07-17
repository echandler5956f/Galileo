#ifndef __galileo_core_constraints_equality_constraint_manager_hpp__
#define __galileo_core_constraints_equality_constraint_manager_hpp__

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "galileo/core/constraints/equality/constraint-generic.hpp"
#include "galileo/core/constraints/equality/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintItemTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ConstraintItemTpl() {}
        ConstraintItemTpl(const std::string &name_, const Model_t &model_, bool active_ = true)
            : name(name_), model(model_), active(active_) {}

        std::string name;
        Model_t model;
        bool active;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = ConstraintCollectionTpl<PS>;
        using ModelManager_t = ConstraintModelManagerTpl<PS, ConstraintCollectionTpl>;
        using DataManager_t = ConstraintDataManagerTpl<PS, ConstraintCollectionTpl>;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Item_t = ConstraintItemTpl<PS, ConstraintCollectionTpl>;

        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;

        using DimNH_t = DimensionTpl<Eigen::Dynamic>;
        static constexpr int NH = DimNH_t::Value;

        using H_t = Eigen::GMatrix<typename PS::VarScalar, NH, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNDX_t::Value, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintDataManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintModelManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintDataManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        ConstraintDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : H(model_manager.get_nh()),
              Hx(model_manager.get_nh(), model_manager.get_ps().ndx_dim.value()),
              Hu(model_manager.get_nh(), model_manager.get_ps().nu_dim.value())
        {
            H.setZero();
            Hx.setZero();
            Hu.setZero();
            for (typename ModelManager_t::ModelContainer_t::const_iterator
                     it = model_manager.get_constraints().begin();
                 it != model_manager.get_constraints().end(); ++it)
            {
                const Item_t &item = it->second;
                constraints.insert(std::make_pair(item.name, item.model.createData(collector)));
            }
        }

        DataContainer_t constraints;
        H_t H;
        Hx_t Hx;
        Hu_t Hu;

    }; // class ConstraintDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintModelManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        ConstraintModelManagerTpl(const PS &ps)
            : ps_(ps), nh_dim_(0)
        {
        }

        void addConstraint(const std::string &name, const Model_t &model, const bool active = true)
        {
            std::pair<typename ModelContainer_t::iterator, bool> ret =
                constraints_.insert(std::make_pair(
                    name, Item_t(name, model, active)));
            if (ret.second == false)
            {
                std::cout << "Warning: we couldn't add the " << name
                          << " constraint item, it already existed." << std::endl;
            }
            else if (active)
            {
                nh_dim_ += model.get_nh_dim();
                active_set_.insert(name);
            }
            else if (!active)
            {
                inactive_set_.insert(name);
            }
        }

        void removeConstraint(const std::string &name)
        {
            typename ModelContainer_t::iterator it = constraints_.find(name);
            if (it != constraints_.end())
            {
                nh_dim_ -= it->second.model.get_nh_dim();
                constraints_.erase(it);
                inactive_set_.erase(name);
            }
            else
            {
                std::cout << "Warning: we couldn't remove the " << name
                          << " constraint item, it doesn't exist." << std::endl;
            }
        }

        void changeConstraintStatus(const std::string &name, bool active)
        {
            typename ModelContainer_t::iterator it = constraints_.find(name);
            if (it != constraints_.end())
            {
                if (active && !it->second.active)
                {
                    nh_dim_ += it->second.model.get_nh_dim();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    nh_dim_ -= it->second.model.get_nh_dim();
                    active_set_.erase(name);
                    inactive_set_.insert(name);
                    it->second.active = active;
                }
            }
            else
            {
                std::cout << "Warning: we couldn't change the status of the " << name
                          << " constraint item, it doesn't exist." << std::endl;
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u)
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived(), u.derived());
                    auto nh_dim_i = m_i.model.get_nh_dim();
                    segment(data.H, nh_accum_i, nh_dim_i) = d_i.H();
                    nh_accum_i += nh_dim_i;
                }
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x)
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived());
                    auto nh_dim_i = m_i.model.get_nh_dim();
                    segment(data.H, nh_accum_i, nh_dim_i) = d_i.H();
                    nh_accum_i += nh_dim_i;
                }
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u)
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived(), u.derived());
                    auto nh_dim_i = m_i.model.get_nh_dim();
                    block(data.Hx, nh_accum_i, 0, nh_dim_i, PS::NDX) = d_i.Hx();
                    block(data.Hu, nh_accum_i, 0, nh_dim_i, PS::NU) = d_i.Hu();
                    nh_accum_i += nh_dim_i;
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived());
                    auto nh_dim_i = m_i.model.get_nh_dim();
                    block(data.Hx, nh_accum_i, 0, nh_dim_i, ps_.ndx_dim) = d_i.Hx();
                    nh_accum_i += nh_dim_i;
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

        bool getConstraintStatus(const std::string &name) const
        {
            typename ModelContainer_t::const_iterator it =
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

        const DimensionTpl<Eigen::Dynamic> &get_nh_dim() const
        {
            return nh_dim_;
        }

        const int get_nh() const
        {
            return nh_dim_.value();
        }

    protected:
        const PS &ps_;
        ModelContainer_t constraints_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

        DimensionTpl<Eigen::Dynamic> nh_dim_;

    }; // class ConstraintModelManagerTpl

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_manager_hpp__
