#ifndef __galileo_multibody_impulses_impulse_manager_hpp__
#define __galileo_multibody_impulses_impulse_manager_hpp__

#include "galileo/common/container/manager-base.hpp"

#include "galileo/multibody/impulses/impulse-base.hpp"
#include "galileo/multibody/impulses/impulse-generic.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseItemTpl
        : public ManagerItemTpl<ImpulseItemTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerItemTpl<ImpulseItemTpl<PhaseSpec, ImpulseCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ImpulseItemTpl(const std::string &name_, const Model_t &model_, const bool active_ = true)
            : Base(name_, model_, active_)
        {
        }

        using Base::active;
        using Base::model;
        using Base::name;

    }; // struct ImpulseItemTpl

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct traits<ImpulseItemTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;
        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct traits<ImpulseManagerTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = ImpulseCollectionTpl<PS>;
        using ModelManager_t = ImpulseModelManagerTpl<PS, ImpulseCollectionTpl>;
        using DataManager_t = ImpulseDataManagerTpl<PS, ImpulseCollectionTpl>;

        using Meta_t = ImpulseTpl<PS, ImpulseCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Item_t = ImpulseItemTpl<PS, ImpulseCollectionTpl>;

        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;

        using Jc_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNV_t::Value, PS::Options>;
        using dv0_dq_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNv_t::Value, PS::Options>;
        using vnext_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, 1, PS::Options>;
        using dnext_dx_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, PS::DimNDX_t::Value, PS::Options>;

        using Force_t = typename PS::Force_t;
        using ForceVector_t = GALILEO_ALIGNED_STD_VECTOR(Force_t);
    };

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct traits<ImpulseDataManagerTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    struct traits<ImpulseModelManagerTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    class ImpulseDataManagerTpl
        : public ManagerDataBase<ImpulseDataManagerTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerDataBase<ImpulseDataManagerTpl<PhaseSpec, ImpulseCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using Jc_t = typename traits<MetaManager_t>::Jc_t;
        using dv0_dq_t = typename traits<MetaManager_t>::dv0_dq_t;
        using vnext_t = typename traits<MetaManager_t>::vnext_t;
        using dnext_dx_t = typename traits<MetaManager_t>::dnext_dx_t;
        using Force_t = typename traits<MetaManager_t>::Force_t;
        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;

        using RobotData_t = typename PS::RobotData_t;

        ImpulseDataManagerTpl(const ModelManager_t &model_manager, RobotData_t *const robot)
            : Base(model_manager, robot),
              fext(model_manager.get_state().get_robot().njoints, Force_t::Zero()),
              Jc(model_manager.get_n_total(), model_manager.get_ps().get_nv()),
              dv0_dq(model_manager.get_n_total(), model_manager.get_ps().get_nv()),
              vnext(model_manager.get_ps().get_nv()),
              dnext_dx(model_manager.get_ps().get_nv(), model_manager.get_ps().get_ndx())
        {
            Jc.setZero();
            dv0_dq.setZero();
            vnext.setZero();
            dnext_dx.setZero();
        }

        using Base::items;
        ForceVector_t fext;
        Jc_t Jc;
        dv0_dq_t dv0_dq;
        vnext_t vnext;
        dnext_dx_t dnext_dx;

    }; // class ImpulseDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    class ImpulseModelManagerTpl
        : public ManagerModelBase<ImpulseModelManagerTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerModelBase<ImpulseModelManagerTpl<PhaseSpec, ImpulseCollectionTpl>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;
        using ForceIterator_t = typename ForceVector_t::iterator;

        ImpulseModelManagerTpl(const PS &ps)
            : Base(),
              ps_(ps), state_(ps.get_state())
        {
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                auto nc_i = get_model_n(m_i.model);
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x);
                    block(data.Jc, nc_accum_i, 0, nc_i, get_ps().get_nv_dim()) = d_i.Jc();
                    nc_accum_i += nc_i;
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                auto nc_i = get_model_n(m_i.model);
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x);
                    block(data.dv0_dq, nc_accum_i, 0, nc_i, get_ps().get_nv_dim()) = d_i.dv0_dq();
                    nc_accum_i += nc_i;
                }
            }
        }

        template <typename VectorNvType>
        void updateVelocity(DataManager_t &data, const Eigen::MatrixBase<VectorNvType> &vnext) const
        {
            data.vnext = vnext;
        }

        template <typename ForceVectorType>
        void updateForce(DataManager_t &data, const Eigen::MatrixBase<ForceVectorType> &force) const
        {
            for (ForceIterator_t it = data.fext.begin(); it != data.fext.end(); ++it)
            {
                *it = PS::Force_t::Zero();
            }

            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                Data_t &d_i = it_d->second;
                auto nc_i = get_model_n(m_i.model);
                if (m_i.active)
                {
                    const auto force_i = segment(force, nc_accum_i, nc_i);
                    m_i.model.updateForce(d_i, force_i);
                    const pinocchio::JointIndex joint =
                        get_state().get_robot().frames[d_i.frame()].parentJoint;
                    data.fext[joint] = d_i.fext();
                    nc_accum_i += nc_i;
                }
                else
                {
                    m_i.model.setZeroForce(d_i);
                }
            }
        }

        template <typename MatrixNvNdxType>
        void updateVelocityDiff(DataManager_t &data, const Eigen::MatrixBase<MatrixNvNdxType> &dnext_dx) const
        {
            data.dnext_dx = dnext_dx;
        }

        template <typename MatrixNcNdxType>
        void updateForceDiff(DataManager_t &data, const Eigen::MatrixBase<MatrixNcNdxType> &df_dx) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                Data_t &d_i = it_d->second;
                auto nc_i = get_model_n(m_i.model);
                if (m_i.active)
                {
                    const auto df_dx_i = block(df_dx, nc_accum_i, 0, nc_i, get_ps().get_ndx_dim());
                    m_i.model.updateForceDiff(d_i, df_dx_i);
                    nc_accum_i += nc_i;
                }
                else
                {
                    m_i.model.setZeroForceDiff(d_i);
                }
            }
        }

        void updateRneaDiff(DataManager_t &data, RobotData_t &robot_data) const
        {
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = items_.begin(), end_m = items_.end(),
                it_d = data.items.begin(), end_d = data.items.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                const Data_t &d_i = it_d->second;
                if (m_i.active)
                {
                    switch (m_i.model.get_type())
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

        template <typename RobotDataType>
        DataManager_t createData(RobotDataType *const robot) const
        {
            return DataManager_t(*this, robot);
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const State_t &get_state() const
        {
            return state_.get();
        }

        int get_model_n(const Model_t &model) const
        {
            return model.get_nc();
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
        std::reference_wrapper<const State_t> state_;

    }; // class ImpulseModelManagerTpl

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_manager_hpp__
