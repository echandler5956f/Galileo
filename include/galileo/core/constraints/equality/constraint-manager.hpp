#ifndef __galileo_core_constraints_equality_constraint_manager_hpp__
#define __galileo_core_constraints_equality_constraint_manager_hpp__

#include "galileo/common/container/manager-base.hpp"

#include "galileo/core/constraints/equality/constraint-base.hpp"
#include "galileo/core/constraints/equality/constraint-generic.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintItemTpl
        : public ManagerItemTpl<ConstraintItemTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerItemTpl<ConstraintItemTpl<PhaseSpec, ConstraintCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ConstraintItemTpl(const std::string &name_, const Model_t &model_, const bool active_ = true)
            : Base(name_, model_, active_)
        {
        }

        using Base::active;
        using Base::model;
        using Base::name;

    }; // struct ConstraintItemTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintItemTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;
        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
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

        using H_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, DimNH_t::Value, PS::DimNU_t::Value, PS::Options>;
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
        : public ManagerDataBase<ConstraintDataManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerDataBase<ConstraintDataManagerTpl<PhaseSpec, ConstraintCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        ConstraintDataManagerTpl(const ModelManager_t &model_manager, DataCollector *const collector)
            : Base(model_manager, collector),
              H(model_manager.get_n_total()),
              Hx(model_manager.get_n_total(), model_manager.get_ps().get_ndx()),
              Hu(model_manager.get_n_total(), model_manager.get_ps().get_nu())
        {
            H.setZero();
            Hx.setZero();
            Hu.setZero();
        }

        using Base::items;
        H_t H;
        Hx_t Hx;
        Hu_t Hu;

    }; // class ConstraintDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintModelManagerTpl
        : public ManagerModelBase<ConstraintModelManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerModelBase<ConstraintModelManagerTpl<PhaseSpec, ConstraintCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        ConstraintModelManagerTpl(const PS &ps)
            : Base(),
              ps_(ps)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x, u);
                    auto nh_i = get_model_n(m_i.model);
                    segment(data.H, nh_accum_i, nh_i) = d_i.H();
                    nh_accum_i += nh_i;
                }
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x);
                    auto nh_i = get_model_n(m_i.model);
                    segment(data.H, nh_accum_i, nh_i) = d_i.H();
                    nh_accum_i += nh_i;
                }
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x, u);
                    auto nh_i = get_model_n(m_i.model);
                    block(data.Hx, nh_accum_i, 0, nh_i, get_ps().get_ndx_dim()) = d_i.Hx();
                    block(data.Hu, nh_accum_i, 0, nh_i, get_ps().get_nu_dim()) = d_i.Hu();
                    nh_accum_i += nh_i;
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nh_accum_i(0);

            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x);
                    auto nh_i = get_model_n(m_i.model);
                    block(data.Hx, nh_accum_i, 0, nh_i, get_ps().get_ndx_dim()) = d_i.Hx();
                    nh_accum_i += nh_i;
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
            return ps_.get();
        }

        int get_model_n(const Model_t &model) const
        {
            return model.get_nh();
        }

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

    }; // class ConstraintModelManagerTpl

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_manager_hpp__
