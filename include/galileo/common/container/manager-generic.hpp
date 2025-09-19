#ifndef __galileo_common_container_manager_generic_hpp__
#define __galileo_common_container_manager_generic_hpp__

#include "galileo/common/container/manager-base.hpp"
#include "galileo/common/container/manager-policies.hpp"

namespace galileo
{
    // =============================================================================
    // GENERIC ITEM TEMPLATE
    // =============================================================================
    
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct GenericItemTpl
        : public ManagerPolicyTraits<Policy>::ItemParams::template ItemImpl<
            ManagerItemTpl<GenericItemTpl<PhaseSpec, CollectionTpl, Policy>>, PhaseSpec>
    {
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Base = typename ManagerPolicyTraits<Policy>::ItemParams::template ItemImpl<
            ManagerItemTpl<GenericItemTpl<PhaseSpec, CollectionTpl, Policy>>, PhaseSpec>;
        
        using Base::Base; // Inherit all constructors
    };

    // =============================================================================
    // GENERIC DATA MANAGER TEMPLATE  
    // =============================================================================
    
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    class GenericDataManagerTpl
        : public ManagerDataBase<GenericDataManagerTpl<PhaseSpec, CollectionTpl, Policy>>
    {
    public:
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerDataBase<GenericDataManagerTpl<PhaseSpec, CollectionTpl, Policy>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;
        using Item_t = typename traits<MetaManager_t>::Item_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        // Constructor - delegates to policy-specific initialization
        template <typename DataCollector>
        GenericDataManagerTpl(const ModelManager_t& model_manager, DataCollector* const collector)
            : Base(model_manager, collector)
        {
            initializeDataMembers(model_manager);
        }

        using Base::items;

        // Policy-specific data members are injected here via template specialization
        // This will be specialized for each concrete policy type

    private:
        template <typename ModelManager>
        void initializeDataMembers(const ModelManager& mm)
        {
            // Default implementation - specializations will provide specific initialization
        }
    };

    // =============================================================================
    // GENERIC MODEL MANAGER TEMPLATE
    // =============================================================================
    
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    class GenericModelManagerTpl
        : public ManagerModelBase<GenericModelManagerTpl<PhaseSpec, CollectionTpl, Policy>>
    {
    public:
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
        using Base = ManagerModelBase<GenericModelManagerTpl<PhaseSpec, CollectionTpl, Policy>>;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;
        using Item_t = typename traits<MetaManager_t>::Item_t;
        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GenericModelManagerTpl(const PS& ps)
            : Base(), ps_(ps)
        {
            initializeManagerSpecific();
        }

        // Calculation methods - delegate to policy
        template <typename StateVectorType, typename... Args>
        void calc(DataManager_t& data, const Eigen::MatrixBase<StateVectorType>& x, Args&&... args) const
        {
            using CalcPolicy = typename ManagerPolicyTraits<Policy>::CalcBehavior;
            if constexpr (requires { CalcPolicy::calc(data, items_, *this, x, std::forward<Args>(args)...); }) {
                CalcPolicy::calc(data, items_, *this, x, std::forward<Args>(args)...);
            } else {
                CalcPolicy::calc(data, items_, x, std::forward<Args>(args)...);
            }
        }

        template <typename StateVectorType, typename... Args>
        void calcDiff(DataManager_t& data, const Eigen::MatrixBase<StateVectorType>& x, Args&&... args) const
        {
            using CalcPolicy = typename ManagerPolicyTraits<Policy>::CalcBehavior;
            if constexpr (requires { CalcPolicy::calcDiff(data, items_, *this, x, std::forward<Args>(args)...); }) {
                CalcPolicy::calcDiff(data, items_, *this, x, std::forward<Args>(args)...);
            } else {
                CalcPolicy::calcDiff(data, items_, x, std::forward<Args>(args)...);
            }
        }

        template <typename DataCollector>
        DataManager_t createData(DataCollector* const collector) const
        {
            return DataManager_t(*this, collector);
        }

        const PS& get_ps() const { return ps_.get(); }

        int get_model_n(const Model_t& model) const
        {
            using ModelDimPolicy = typename ManagerPolicyTraits<Policy>::ModelDimension;
            return ModelDimPolicy::get_n(model);
        }

        // Inherit standard interface from base
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

        // Hook for policy-specific initialization
        void initializeManagerSpecific()
        {
            // Default implementation - can be specialized for specific policies
        }
    };

    // =============================================================================
    // GENERIC MANAGER META-TEMPLATE
    // =============================================================================
    
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct GenericManagerTpl
    {
        using PS = PhaseSpec;
        using Policy_t = Policy;
    };

    // =============================================================================
    // AUTOMATIC TRAITS SPECIALIZATION GENERATOR
    // =============================================================================
    
    /**
     * @brief Macro-free way to generate all required traits specializations
     * Just specialize this template for each concrete policy
     */
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct traits<GenericItemTpl<PhaseSpec, CollectionTpl, Policy>>
    {
        using MetaManager_t = GenericManagerTpl<PhaseSpec, CollectionTpl, Policy>;
    };

    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct traits<GenericManagerTpl<PhaseSpec, CollectionTpl, Policy>>
        : public GenerateManagerTraits<PhaseSpec, CollectionTpl, Policy,
                                     GenericItemTpl, GenericDataManagerTpl, GenericModelManagerTpl> {};

    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct traits<GenericDataManagerTpl<PhaseSpec, CollectionTpl, Policy>>
    {
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct traits<GenericModelManagerTpl<PhaseSpec, CollectionTpl, Policy>>
    {
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

} // namespace galileo

#endif // __galileo_common_container_manager_generic_hpp__