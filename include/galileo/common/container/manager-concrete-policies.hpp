#ifndef __galileo_common_container_manager_concrete_policies_hpp__
#define __galileo_common_container_manager_concrete_policies_hpp__

#include "galileo/common/container/manager-generic.hpp"

namespace galileo
{
    // =============================================================================
    // COST MANAGER POLICY
    // =============================================================================
    
    template <typename PhaseSpec>
    struct CostManagerPolicy
    {
        using PS = PhaseSpec;
        using NumScalar = typename PS::NumScalar;
        
        using ItemParams = WeightedItemParams<NumScalar>;
        using DataMembers = void; // Will be specialized below
        using CalcBehavior = AccumulationCalc;
        using ModelDimension = ZeroDimension;
    };

    // Specialized data manager for cost policy
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    class GenericDataManagerTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>
        : public ManagerDataBase<GenericDataManagerTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>>
    {
    public:
        using PS = PhaseSpec;
        using Policy = CostManagerPolicy<PS>;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Base = ManagerDataBase<GenericDataManagerTpl<PS, CollectionTpl, Policy>>;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;

        // Cost-specific typedefs (maintaining compatibility)
        GALILEO_COST_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        GenericDataManagerTpl(const ModelManager_t& model_manager, DataCollector* const collector)
            : Base(model_manager, collector),
              L(0.),
              Lx(model_manager.get_ps().get_ndx()),
              Lu(model_manager.get_ps().get_nu()),
              Lxx(model_manager.get_ps().get_ndx(), model_manager.get_ps().get_ndx()),
              Lxu(model_manager.get_ps().get_ndx(), model_manager.get_ps().get_nu()),
              Luu(model_manager.get_ps().get_nu(), model_manager.get_ps().get_nu())
        {
            Lx.setZero(); Lu.setZero(); Lxx.setZero(); Lxu.setZero(); Luu.setZero();
        }

        using Base::items;
        L_t L;
        Lx_t Lx;
        Lu_t Lu;
        Lxx_t Lxx;
        Lxu_t Lxu;
        Luu_t Luu;
    };

    // =============================================================================
    // CONSTRAINT MANAGER POLICY
    // =============================================================================
    
    template <typename PhaseSpec>
    struct ConstraintManagerPolicy
    {
        using PS = PhaseSpec;
        
        using ItemParams = StandardItemParams;
        using DataMembers = void; // Will be specialized below
        using CalcBehavior = SegmentationCalc;
        using ModelDimension = ConstraintDimension;
    };

    // Specialized data manager for constraint policy
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    class GenericDataManagerTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>
        : public ManagerDataBase<GenericDataManagerTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>>
    {
    public:
        using PS = PhaseSpec;
        using Policy = ConstraintManagerPolicy<PS>;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Base = ManagerDataBase<GenericDataManagerTpl<PS, CollectionTpl, Policy>>;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;

        // Constraint-specific typedefs (maintaining compatibility)
        GALILEO_CONSTRAINT_DATA_TYPEDEF(MetaManager_t);

        template <typename DataCollector>
        GenericDataManagerTpl(const ModelManager_t& model_manager, DataCollector* const collector)
            : Base(model_manager, collector),
              H(model_manager.get_n_total()),
              Hx(model_manager.get_n_total(), model_manager.get_ps().get_ndx()),
              Hu(model_manager.get_n_total(), model_manager.get_ps().get_nu())
        {
            H.setZero(); Hx.setZero(); Hu.setZero();
        }

        using Base::items;
        H_t H;
        Hx_t Hx;
        Hu_t Hu;
    };

    // =============================================================================
    // CONTACT MANAGER POLICY
    // =============================================================================
    
    template <typename PhaseSpec>
    struct ContactManagerPolicy
    {
        using PS = PhaseSpec;
        
        using ItemParams = StandardItemParams;
        using DataMembers = void; // Will be specialized below  
        using CalcBehavior = ContactCalc; // Will define this below
        using ModelDimension = ContactDimension;
    };

    // Contact-specific calculation behavior
    struct ContactCalc
    {
        template <typename DataManager, typename ModelContainer, typename ModelManager,
                  typename StateVectorType>
        static void calc(DataManager& data, const ModelContainer& items, const ModelManager& mm,
                        const Eigen::MatrixBase<StateVectorType>& x)
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            
            for (auto it_m = items.begin(), it_d = data.items.begin();
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calc(d_i, x);
                    auto nc_i = mm.get_model_n(m_i.model);
                    segment(data.a0, nc_accum_i, nc_i) = d_i.a0();
                    block(data.Jc, nc_accum_i, 0, nc_i, mm.get_ps().get_nv_dim()) = d_i.Jc();
                    nc_accum_i += nc_i;
                }
            }
        }

        template <typename DataManager, typename ModelContainer, typename ModelManager,
                  typename StateVectorType>  
        static void calcDiff(DataManager& data, const ModelContainer& items, const ModelManager& mm,
                           const Eigen::MatrixBase<StateVectorType>& x)
        {
            DimensionTpl<Eigen::Dynamic> nc_accum_i(0);
            
            for (auto it_m = items.begin(), it_d = data.items.begin();
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calcDiff(d_i, x);
                    auto nc_i = mm.get_model_n(m_i.model);
                    block(data.da0_dx, nc_accum_i, 0, nc_i, mm.get_ps().get_ndx_dim()) = d_i.da0_dx();
                    nc_accum_i += nc_i;
                }
            }
        }
    };

    // Specialized data manager for contact policy
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    class GenericDataManagerTpl<PhaseSpec, CollectionTpl, ContactManagerPolicy<PhaseSpec>>
        : public ManagerDataBase<GenericDataManagerTpl<PhaseSpec, CollectionTpl, ContactManagerPolicy<PhaseSpec>>>
    {
    public:
        using PS = PhaseSpec;
        using Policy = ContactManagerPolicy<PS>;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Base = ManagerDataBase<GenericDataManagerTpl<PS, CollectionTpl, Policy>>;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using RobotData_t = typename PS::RobotData_t;

        // Contact-specific types
        using Force_t = typename PS::Force_t;
        using ForceVector_t = GALILEO_ALIGNED_STD_VECTOR(Force_t);
        using Jc_t = typename traits<MetaManager_t>::Jc_t;
        using a0_t = typename traits<MetaManager_t>::a0_t;
        using da0_dx_t = typename traits<MetaManager_t>::da0_dx_t;
        using dv_t = typename traits<MetaManager_t>::dv_t;
        using ddv_dx_t = typename traits<MetaManager_t>::ddv_dx_t;

        GenericDataManagerTpl(const ModelManager_t& model_manager, RobotData_t* const robot)
            : Base(model_manager, robot),
              fext(model_manager.get_state().get_robot().njoints, Force_t::Zero()),
              Jc(model_manager.get_n_total(), model_manager.get_ps().get_nv()),
              a0(model_manager.get_n_total()),
              da0_dx(model_manager.get_n_total(), model_manager.get_ps().get_ndx()),
              dv(model_manager.get_ps().get_nv()),
              ddv_dx(model_manager.get_ps().get_nv(), model_manager.get_ps().get_ndx())
        {
            Jc.setZero(); a0.setZero(); da0_dx.setZero(); dv.setZero(); ddv_dx.setZero();
        }

        using Base::items;
        ForceVector_t fext;
        Jc_t Jc;
        a0_t a0;
        da0_dx_t da0_dx;
        dv_t dv;
        ddv_dx_t ddv_dx;
    };

    // Specialized model manager for contact policy (adds extra methods)
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    class GenericModelManagerTpl<PhaseSpec, CollectionTpl, ContactManagerPolicy<PhaseSpec>>
        : public ManagerModelBase<GenericModelManagerTpl<PhaseSpec, CollectionTpl, ContactManagerPolicy<PhaseSpec>>>
    {
    public:
        using PS = PhaseSpec;
        using Policy = ContactManagerPolicy<PS>;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Base = ManagerModelBase<GenericModelManagerTpl<PS, CollectionTpl, Policy>>;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using RobotData_t = typename PS::RobotData_t;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        GenericModelManagerTpl(const PS& ps)
            : Base(), ps_(ps), state_(ps.get_state()) {}

        // Standard calc/calcDiff from base + policy
        using Base::calc;
        using Base::calcDiff;

        // Contact-specific additional methods
        template <typename VectorNvType>
        void updateAcceleration(DataManager_t& data, const Eigen::MatrixBase<VectorNvType>& dv) const
        {
            data.dv = dv;
        }

        template <typename ForceVectorType>
        void updateForce(DataManager_t& data, const Eigen::MatrixBase<ForceVectorType>& force) const
        {
            // Implementation similar to original ContactModelManagerTpl
            // ... [detailed implementation would go here]
        }

        // ... other contact-specific methods

        const State_t& get_state() const { return state_.get(); }

        using Base::get_ps;
        using Base::get_model_n;
        using Base::addItem; using Base::removeItem; using Base::changeItemStatus;
        using Base::get_active_set; using Base::get_inactive_set; using Base::get_items;
        using Base::get_item_status; using Base::get_n_active; using Base::get_n_active_dim;
        using Base::get_n_total; using Base::get_n_total_dim;

    protected:
        using Base::items_; using Base::active_set_; using Base::inactive_set_;
        using Base::active_dim_; using Base::total_dim_;
        std::reference_wrapper<const PS> ps_;
        std::reference_wrapper<const State_t> state_;
    };

    // =============================================================================
    // IMPULSE MANAGER POLICY (Similar to Contact)
    // =============================================================================
    
    template <typename PhaseSpec>
    struct ImpulseManagerPolicy
    {
        using PS = PhaseSpec;
        
        using ItemParams = StandardItemParams;
        using DataMembers = void; // Will be specialized
        using CalcBehavior = ImpulseCalc; // Similar to ContactCalc
        using ModelDimension = ContactDimension; // Same as contact (uses get_nc())
    };

    // ... [Similar specializations for ImpulseManagerPolicy would follow]

    // =============================================================================
    // BACKWARD COMPATIBILITY TYPE ALIASES
    // =============================================================================
    
    // These maintain full backward compatibility with existing code
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using CostManagerTpl = GenericManagerTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using CostItemTpl = GenericItemTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using CostDataManagerTpl = GenericDataManagerTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using CostModelManagerTpl = GenericModelManagerTpl<PhaseSpec, CollectionTpl, CostManagerPolicy<PhaseSpec>>;

    // Similar aliases for other managers...
    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using ConstraintManagerTpl = GenericManagerTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using ConstraintItemTpl = GenericItemTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using ConstraintDataManagerTpl = GenericDataManagerTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>;

    template <typename PhaseSpec, template <typename PS> class CollectionTpl>
    using ConstraintModelManagerTpl = GenericModelManagerTpl<PhaseSpec, CollectionTpl, ConstraintManagerPolicy<PhaseSpec>>;

    // ... and so on for Contact and Impulse managers

} // namespace galileo

#endif // __galileo_common_container_manager_concrete_policies_hpp__