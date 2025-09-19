#ifndef __galileo_common_container_manager_policies_hpp__
#define __galileo_common_container_manager_policies_hpp__

#include "galileo/fwd.hpp"
#include <type_traits>

namespace galileo
{
    // Forward declarations for policy-based manager generation
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy>
    struct GenericManagerTpl;

    // =============================================================================
    // MANAGER POLICY CONCEPT
    // =============================================================================
    
    /**
     * @brief Policy concept for defining manager behavior variations
     * 
     * A policy must provide:
     * - ItemParams: tuple of additional constructor parameters beyond (name, model, active)
     * - DataMembers: type list defining data manager member variables  
     * - CalcBehavior: calculation implementation policy
     * - ModelDimension: policy for getting model dimensions
     */
    template <typename Policy>
    struct ManagerPolicyTraits
    {
        using ItemParams = typename Policy::ItemParams;
        using DataMembers = typename Policy::DataMembers;
        using CalcBehavior = typename Policy::CalcBehavior;
        using ModelDimension = typename Policy::ModelDimension;
    };

    // =============================================================================
    // ITEM PARAMETER POLICIES
    // =============================================================================
    
    // No additional parameters beyond standard (name, model, active)
    struct StandardItemParams 
    {
        using type = std::tuple<>;
        
        template <typename Base, typename PS>
        struct ItemImpl : public Base
        {
            using Base::Base; // Inherit constructors
        };
    };

    // Add weight parameter for cost items
    template <typename NumScalar>
    struct WeightedItemParams
    {
        using type = std::tuple<NumScalar>;
        
        template <typename Base, typename PS>
        struct ItemImpl : public Base
        {
            using typename Base::Model_t;
            using NumScalar_t = NumScalar;
            
            ItemImpl(const std::string& name_, const Model_t& model_, const NumScalar_t& weight_, bool active_ = true)
                : Base(name_, model_, active_), weight(weight_) {}
                
            using Base::active;
            using Base::model; 
            using Base::name;
            NumScalar_t weight;
        };
    };

    // =============================================================================
    // DATA MEMBER SPECIFICATION
    // =============================================================================
    
    // Helper for specifying data member types and initialization
    template <typename... Members>
    struct DataMemberSpec
    {
        using type = std::tuple<Members...>;
    };
    
    // Individual data member specification
    template <typename T, typename InitFunc>
    struct DataMember
    {
        using type = T;
        using init_func = InitFunc;
    };

    // Initialization functions for common patterns
    struct ZeroInit 
    {
        template <typename T, typename ModelManager>
        static T create(const ModelManager& mm) 
        {
            auto result = T(mm.get_n_total());
            result.setZero();
            return result;
        }
    };

    struct MatrixInit
    {
        template <typename T, typename ModelManager>
        static T create(const ModelManager& mm, int rows, int cols)
        {
            auto result = T(rows, cols);
            result.setZero();
            return result;
        }
    };

    // =============================================================================
    // CALCULATION BEHAVIOR POLICIES
    // =============================================================================
    
    // Accumulation pattern (used by CostManager)
    struct AccumulationCalc
    {
        template <typename DataManager, typename ModelContainer, typename StateVectorType, typename... Args>
        static void calc(DataManager& data, const ModelContainer& items, 
                        const Eigen::MatrixBase<StateVectorType>& x, Args&&... args)
        {
            using VarScalar = typename DataManager::MetaManager_t::PS::VarScalar;
            data.L = VarScalar(0.);
            
            for (auto it_m = items.begin(), it_d = data.items.begin(); 
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calc(d_i, x, std::forward<Args>(args)...);
                    data.L += m_i.weight * d_i.L();
                }
            }
        }
        
        template <typename DataManager, typename ModelContainer, typename StateVectorType, typename... Args>
        static void calcDiff(DataManager& data, const ModelContainer& items,
                           const Eigen::MatrixBase<StateVectorType>& x, Args&&... args)
        {
            data.Lx.setZero(); data.Lu.setZero(); 
            data.Lxx.setZero(); data.Lxu.setZero(); data.Luu.setZero();
            
            for (auto it_m = items.begin(), it_d = data.items.begin();
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calcDiff(d_i, x, std::forward<Args>(args)...);
                    data.Lx += m_i.weight * d_i.Lx();
                    data.Lu += m_i.weight * d_i.Lu();
                    data.Lxx += m_i.weight * d_i.Lxx(); 
                    data.Lxu += m_i.weight * d_i.Lxu();
                    data.Luu += m_i.weight * d_i.Luu();
                }
            }
        }
    };

    // Segmentation pattern (used by ConstraintManager) 
    struct SegmentationCalc
    {
        template <typename DataManager, typename ModelContainer, typename ModelManager, 
                  typename StateVectorType, typename... Args>
        static void calc(DataManager& data, const ModelContainer& items, const ModelManager& mm,
                        const Eigen::MatrixBase<StateVectorType>& x, Args&&... args)
        {
            DimensionTpl<Eigen::Dynamic> accum_i(0);
            
            for (auto it_m = items.begin(), it_d = data.items.begin();
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calc(d_i, x, std::forward<Args>(args)...);
                    auto n_i = mm.get_model_n(m_i.model);
                    segment(data.H, accum_i, n_i) = d_i.H();
                    accum_i += n_i;
                }
            }
        }
        
        template <typename DataManager, typename ModelContainer, typename ModelManager,
                  typename StateVectorType, typename... Args>  
        static void calcDiff(DataManager& data, const ModelContainer& items, const ModelManager& mm,
                           const Eigen::MatrixBase<StateVectorType>& x, Args&&... args)
        {
            DimensionTpl<Eigen::Dynamic> accum_i(0);
            
            for (auto it_m = items.begin(), it_d = data.items.begin();
                 it_m != items.end(); ++it_m, ++it_d)
            {
                const auto& m_i = it_m->second;
                if (m_i.active)
                {
                    auto& d_i = it_d->second;
                    m_i.model.calcDiff(d_i, x, std::forward<Args>(args)...);
                    auto n_i = mm.get_model_n(m_i.model);
                    block(data.Hx, accum_i, 0, n_i, mm.get_ps().get_ndx_dim()) = d_i.Hx();
                    if constexpr (sizeof...(args) > 0) {
                        block(data.Hu, accum_i, 0, n_i, mm.get_ps().get_nu_dim()) = d_i.Hu();
                    }
                    accum_i += n_i;
                }
            }
        }
    };

    // =============================================================================
    // MODEL DIMENSION POLICIES
    // =============================================================================
    
    struct ZeroDimension
    {
        template <typename Model>
        static int get_n(const Model& model) { return 0; }
    };
    
    struct ConstraintDimension  
    {
        template <typename Model>
        static int get_n(const Model& model) { return model.get_nh(); }
    };
    
    struct ContactDimension
    {
        template <typename Model> 
        static int get_n(const Model& model) { return model.get_nc(); }
    };

    // =============================================================================
    // AUTOMATIC TRAITS GENERATION
    // =============================================================================
    
    /**
     * @brief Generates all required traits for a manager given a policy
     */
    template <typename PhaseSpec, template <typename PS> class CollectionTpl, typename Policy,
              typename ItemTemplate = GenericItemTpl,
              typename DataTemplate = GenericDataManagerTpl,
              typename ModelTemplate = GenericModelManagerTpl>
    struct GenerateManagerTraits
    {
        using PS = PhaseSpec;
        using MetaManager_t = GenericManagerTpl<PS, CollectionTpl, Policy>;
        using Collection_t = CollectionTpl<PS>;
        using ModelManager_t = ModelTemplate<PS, CollectionTpl, Policy>;
        using DataManager_t = DataTemplate<PS, CollectionTpl, Policy>;
        
        using Meta_t = typename CollectionTpl<PS>::template Manager<Policy>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        
        using Item_t = ItemTemplate<PS, CollectionTpl, Policy>;
        
        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;
    };

} // namespace galileo

#endif // __galileo_common_container_manager_policies_hpp__