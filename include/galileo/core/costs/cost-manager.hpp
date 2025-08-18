#ifndef __galileo_core_costs_cost_manager_hpp__
#define __galileo_core_costs_cost_manager_hpp__

#include "galileo/common/container/manager-base.hpp"

#include "galileo/core/costs/cost-base.hpp"
#include "galileo/core/costs/cost-generic.hpp"

// Despite its name, this file only pertains to cost residuals,
// not just any general costs derived from CostModelBase/CostDataBase

namespace galileo
{

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostManagerTpl;

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostItemTpl : public ManagerItemTpl<CostItemTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerItemTpl<CostItemTpl<PhaseSpec, CostCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using NumScalar = typename PS::NumScalar;

        CostItemTpl(const std::string &name_,
                    const Model_t &model_,
                    const NumScalar &weight_,
                    const bool active_ = true)
            : Base(name_, model_, active_), weight(weight_)
        {
        }

        using Base::active;
        using Base::model;
        using Base::name;
        NumScalar weight;

    }; // struct CostItemTpl

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostItemTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;
        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
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

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostDataManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct traits<CostModelManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    class CostDataManagerTpl : public ManagerDataBase<CostDataManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerDataBase<CostDataManagerTpl<PhaseSpec, CostCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GALILEO_COST_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        CostDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : Base(model_manager, collector),
              L(0.),
              Lx(model_manager.get_ps().get_ndx()),
              Lu(model_manager.get_ps().get_nu()),
              Lxx(model_manager.get_ps().get_ndx(), model_manager.get_ps().get_ndx()),
              Lxu(model_manager.get_ps().get_ndx(), model_manager.get_ps().get_nu()),
              Luu(model_manager.get_ps().get_nu(), model_manager.get_ps().get_nu())
        {
            Lx.setZero();
            Lu.setZero();
            Lxx.setZero();
            Lxu.setZero();
            Luu.setZero();
        }

        using Base::items;
        L_t L;
        Lx_t Lx;
        Lu_t Lu;
        Lxx_t Lxx;
        Lxu_t Lxu;
        Luu_t Luu;

    }; // class CostDataManagerTpl

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    class CostModelManagerTpl : public ManagerModelBase<CostModelManagerTpl<PhaseSpec, CostCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = CostManagerTpl<PS, CostCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerModelBase<CostModelManagerTpl<PhaseSpec, CostCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using NumScalar = typename PS::NumScalar;
        using VarScalar = typename PS::VarScalar;

        CostModelManagerTpl(const PS &ps) : Base(), ps_(ps) {}

        template <typename StateVectorType, typename ControlVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            data.L = VarScalar(0.);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(), it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d;
                 ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x, u);
                    data.L += m_i.weight * d_i.L();
                }
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            data.L = VarScalar(0.);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(), it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d;
                 ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x);
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

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(), it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d;
                 ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x, u);
                    data.Lx += m_i.weight * d_i.Lx();
                    data.Lu += m_i.weight * d_i.Lu();
                    data.Lxx += m_i.weight * d_i.Lxx();
                    data.Lxu += m_i.weight * d_i.Lxu();
                    data.Luu += m_i.weight * d_i.Luu();
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            data.Lx.setZero();
            data.Lxx.setZero();

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(), it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d;
                 ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x);
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

        const PS &get_ps() const { return ps_.get(); }
        int get_model_n(const Model_t &model) const { return 0; }

        using Base::addItem;
        using Base::removeItem;
        using Base::changeItemStatus;
        using Base::get_active_set;
        using Base::get_inactive_set;
        using Base::get_items;
        using Base::get_item_status;
        using Base::get_n_active;
        using Base::get_n_active_dim;
        using Base::get_n_total;
        using Base::get_n_total_dim;

    protected:
        using Base::items_;
        using Base::active_set_;
        using Base::inactive_set_;
        using Base::active_dim_;
        using Base::total_dim_;

        std::reference_wrapper<const PS> ps_;

    }; // class CostModelManagerTpl

} // namespace galileo

#endif // __galileo_core_costs_cost_manager_hpp__
