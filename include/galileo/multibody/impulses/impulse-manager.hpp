#ifndef __galileo_multibody_impulses_impulse_manager_hpp__
#define __galileo_multibody_impulses_impulse_manager_hpp__

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "galileo/multibody/impulses/fwd.hpp"

#include "galileo/multibody/force-base.hpp"
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
    {
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ImpulseItemTpl() {}
        ImpulseItemTpl(const std::string &name_, const Model_t &model_, const bool active_ = true)
            : name(name_), model(model_), active(active_) {}

        std::string name;
        Model_t model;
        bool active;
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
    {
    public:
        using PS = PhaseSpec;

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

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
            : fext(model_manager.get_state().get_robot().njoints, Force_t::Zero()),
              Jc(model_manager.get_nc_total(), model_manager.get_ps().get_nv()),
              dv0_dq(model_manager.get_nc_total(), model_manager.get_ps().get_nv()),
              vnext(model_manager.get_ps().get_nv()),
              dnext_dx(model_manager.get_ps().get_nv(), model_manager.get_ps().get_ndx())
        {
            Jc.setZero();
            dv0_dq.setZero();
            vnext.setZero();
            dnext_dx.setZero();
            for (typename ModelContainer_t::const_iterator
                     it = model_manager.get_impulses().begin();
                 it != model_manager.get_impulses().end(); ++it)
            {
                const Item_t &item = it->second;
                impulses.insert(
                    std::make_pair(item.name, item.model.createData(robot)));
            }
        }

        DataContainer_t impulses;
        ForceVector_t fext;

        Jc_t Jc;
        dv0_dq_t dv0_dq;
        vnext_t vnext;
        dnext_dx_t dnext_dx;

    }; // class ImpulseDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ImpulseCollectionTpl>
    class ImpulseModelManagerTpl
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using MetaManager_t = ImpulseManagerTpl<PS, ImpulseCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using ForceVector_t = typename traits<MetaManager_t>::ForceVector_t;
        using ForceIterator_t = typename ForceVector_t::iterator;

        ImpulseModelManagerTpl(const PS &ps)
            : ps_(ps), state_(ps.get_state()),
              nc_active_dim_(DimensionTpl<Eigen::Dynamic>(0)), nc_total_dim_(DimensionTpl<Eigen::Dynamic>(0))
        {
        }

        void addImpulse(const std::string &name, const Model_t &model, bool active = true)
        {
            std::pair<typename ModelContainer_t::iterator, bool> ret =
                impulses_.insert(std::make_pair(
                    name, Item_t(name, model, active)));
            if (ret.second == false)
            {
                std::cerr << "Warning: we couldn't add the " << name
                          << " impulse item, it already existed." << std::endl;
            }
            else if (active)
            {
                nc_active_dim_ += model.get_nc();
                nc_total_dim_ += model.get_nc();
                active_set_.insert(name);
            }
            else if (!active)
            {
                nc_total_dim_ += model.get_nc();
                inactive_set_.insert(name);
            }
        }

        void removeImpulse(const std::string &name)
        {
            typename ModelContainer_t::iterator it = impulses_.find(name);
            if (it != impulses_.end())
            {
                nc_active_dim_ -= it->second.model.get_nc();
                nc_total_dim_ -= it->second.model.get_nc();
                impulses_.erase(it);
                inactive_set_.erase(name);
            }
            else
            {
                std::cerr << "Warning: we couldn't remove the " << name
                          << " impulse item, it doesn't exist." << std::endl;
            }
        }

        void changeImpulseStatus(const std::string &name, bool active)
        {
            typename ModelContainer_t::iterator it = impulses_.find(name);
            if (it != impulses_.end())
            {
                if (active && !it->second.active)
                {
                    nc_active_dim_ += it->second.model.get_nc();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    nc_active_dim_ -= it->second.model.get_nc();
                    active_set_.erase(name);
                    inactive_set_.insert(name);
                    it->second.active = active;
                }
                it->second.active = active;
            }
            else
            {
                std::cerr << "Warning: we couldn't change the status of the " << name
                          << " impulse item, it doesn't exist." << std::endl;
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = impulses_.begin(), end_m = impulses_.end(),
                it_d = data.impulses.begin(), end_d = data.impulses.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                auto nc_dim_i = m_i.model.get_nc();
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x);
                    block(data.Jc, nc_accum_i, 0, nc_dim_i, get_ps().get_nv_dim()) = d_i.Jc();
                    nc_accum_i += nc_dim_i;
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            typename ModelContainer_t::const_iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = impulses_.begin(), end_m = impulses_.end(),
                it_d = data.impulses.begin(), end_d = data.impulses.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                auto nc_dim_i = m_i.model.get_nc();
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x);
                    block(data.dv0_dq, nc_accum_i, 0, nc_dim_i, get_ps().get_nv_dim()) = d_i.dv0_dq();
                    nc_accum_i += nc_dim_i;
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
            for (it_m = impulses_.begin(), end_m = impulses_.end(),
                it_d = data.impulses.begin(), end_d = data.impulses.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                Data_t &d_i = it_d->second;
                auto nc_dim_i = m_i.model.get_nc();
                if (m_i.active)
                {
                    const auto force_i = segment(force, nc_accum_i, nc_dim_i);
                    m_i.model.updateForce(d_i, force_i);
                    const pinocchio::JointIndex joint =
                        get_state().get_robot().frames[d_i.frame()].parentJoint;
                    data.fext[joint] = d_i.fext();
                    nc_accum_i += nc_dim_i;
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
            for (it_m = impulses_.begin(), end_m = impulses_.end(),
                it_d = data.impulses.begin(), end_d = data.impulses.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                Data_t &d_i = it_d->second;
                auto nc_dim_i = m_i.model.get_nc();
                if (m_i.active)
                {
                    const auto df_dx_i = block(df_dx, nc_accum_i, 0, nc_dim_i, get_ps().get_ndx_dim());
                    m_i.model.updateForceDiff(d_i, df_dx_i);
                    nc_accum_i += nc_dim_i;
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
            for (it_m = impulses_.begin(), end_m = impulses_.end(),
                it_d = data.impulses.begin(), end_d = data.impulses.end();
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

        DataManager_t createData(RobotData_t *const robot) const
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

        const ModelContainer_t &get_impulses() const
        {
            return impulses_;
        }

        const DimensionTpl<Eigen::Dynamic> &get_nc_active_dim() const
        {
            return nc_active_dim_;
        }

        int get_nc_active() const
        {
            return nc_active_dim_.value();
        }

        const DimensionTpl<Eigen::Dynamic> &get_nc_total_dim() const
        {
            return nc_total_dim_;
        }

        int get_nc_total() const
        {
            return nc_total_dim_.value();
        }

    protected:
        std::reference_wrapper<const PS> ps_;
        std::reference_wrapper<const State_t> state_;
        ModelContainer_t impulses_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

        DimensionTpl<Eigen::Dynamic> nc_active_dim_;
        DimensionTpl<Eigen::Dynamic> nc_total_dim_;

    }; // class ImpulseModelManagerTpl

} // namespace galileo

#endif // __galileo_multibody_impulses_impulse_manager_hpp__
