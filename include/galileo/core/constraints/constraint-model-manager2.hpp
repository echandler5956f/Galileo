#ifndef __galileo_core_constraints_constraint_model_manager_hpp__
#define __galileo_core_constraints_constraint_model_manager_hpp__

#include <tuple>
#include <array>
#include <cstdint>
#include <initializer_list>

namespace galileo
{

    namespace core
    {

        //-----------------------------------------------------------------
        // Helper: Compute the index of a type T in a tuple.
        template <typename T, typename Tuple>
        struct tuple_index;

        template <typename T, typename... Ts>
        struct tuple_index<T, std::tuple<T, Ts...>>
        {
            static constexpr std::size_t value = 0;
        }; // struct tuple_index

        template <typename T, typename U, typename... Ts>
        struct tuple_index<T, std::tuple<U, Ts...>>
        {
            static constexpr std::size_t value = 1 + tuple_index<T, std::tuple<Ts...>>::value;
        }; // struct tuple_index

        //-----------------------------------------------------------------
        // ConstraintModelManager:
        // - Holds a tuple of constraint models and an active bitmask.
        // - Provides type–safe activate()/deactivate() functions.
        // - Before calling calcAll(), the user should call updateDataCollectionSizes()
        //   to resize g_active and h_active to exactly fit the sum of active constraints’ sizes.
        // - In calcAll(), if all constraints are active, a fast, fully unrolled code path is used;
        //   otherwise, a cached code path that tests each model’s active flag is used.
        template <typename... ConstraintModels>
        class ConstraintModelManager
        {
        public:
            std::tuple<ConstraintModels...> models_;
            // Bit mask: bit i corresponds to the i-th constraint in the tuple.
            // A bit set to 1 means that constraint is active.
            uint64_t active_mask_;

            // full_mask_ is all ones (for the number of models).
            static constexpr uint64_t full_mask_ = (sizeof...(ConstraintModels) >= 64 ? ~0ULL : ((1ULL << sizeof...(ConstraintModels)) - 1));

            // Constructor: all constraints active by default.
            template <typename... Models>
            ConstraintModelManager(Models &&...models)
                : models_(std::forward<Models>(models)...),
                  active_mask_(full_mask_)
            {
            }

            // --- Single-constraint activation:
            template <typename ConstraintType>
            void activate()
            {
                constexpr std::size_t idx = tuple_index<ConstraintType, std::tuple<ConstraintModels...>>::value;
                active_mask_ |= (1ULL << idx);
            }

            // Variadic activation (only enabled when at least 2 types are provided):
            template <typename First, typename Second, typename... Others>
            void activate()
            {
                activate<First>();
                activate<Second, Others...>();
            }

            // --- Single-constraint deactivation:
            template <typename ConstraintType>
            void deactivate()
            {
                constexpr std::size_t idx = tuple_index<ConstraintType, std::tuple<ConstraintModels...>>::value;
                active_mask_ &= ~(1ULL << idx);
            }

            // Variadic deactivation (only enabled when at least 2 types are provided):
            template <typename First, typename Second, typename... Others>
            void deactivate()
            {
                deactivate<First>();
                deactivate<Second, Others...>();
            }

            // updateDataCollectionSizes:
            // Walks over the models (compile-time unrolled) and sums the sizes of the g and h
            // buffers for all active constraints. Then resizes data_collection accordingly.
            void updateDataCollectionSizes(ConstraintDataCollection &data_collection) const
            {
                int total_g = 0;
                int total_h = 0;
                updateDataCollectionSizesImpl<0>(total_g, total_h);
                data_collection.g_active.resize(total_g);
                data_collection.h_active.resize(total_h);
            }

            // calcAll:
            // Calls each active constraint’s calc() and assigns its output buffers into the
            // corresponding segments of data_collection.
            // Two code paths are provided:
            //  - Fast path: when all constraints are active.
            //  - Cached path: when only a subset is active.
            template <typename StateVector, typename ControlVector>
            void calcAll(ConstraintDataCollection &data_collection,
                         const Eigen::MatrixBase<StateVector> &x,
                         const Eigen::MatrixBase<ControlVector> &u)
            {
                if (active_mask_ == full_mask_)
                {
                    // Fast path: all constraints active.
                    calcAllFastPath(x.derived(), u.derived(), data_collection);
                }
                else
                {
                    // Cached path: only some constraints active.
                    calcAllCachedPath(x.derived(), u.derived(), data_collection);
                }
            }

        private:
            //-----------------------------------------------------------------
            // Helpers to compute compile-time cumulative offsets.
            // For each index I, these return the starting offset (in g or h)
            // for the I-th constraint.
            template <std::size_t I>
            static constexpr int g_offset()
            {
                if constexpr (I == 0)
                    return 0;
                else
                    return g_offset<I - 1>() + std::tuple_element_t<I - 1, std::tuple<ConstraintModels...>>::DataType::G_SIZE;
            }

            template <std::size_t I>
            static constexpr int h_offset()
            {
                if constexpr (I == 0)
                    return 0;
                else
                    return h_offset<I - 1>() + std::tuple_element_t<I - 1, std::tuple<ConstraintModels...>>::DataType::H_SIZE;
            }

            //-----------------------------------------------------------------
            // Fast path implementation (all constraints active):
            // Unrolls the calls for each constraint using an index sequence.
            template <typename StateVector, typename ControlVector, std::size_t... Is>
            void calcAllFastPathImpl(const StateVector &x,
                                     const ControlVector &u,
                                     ConstraintDataCollection &data_collection,
                                     std::index_sequence<Is...>)
            {
                (void)std::initializer_list<int>{
                    ([&]()
                     {
          using ModelType = std::tuple_element_t<Is, std::tuple<ConstraintModels...>>;
          using DataType = typename ModelType::DataType;
          DataType local_data;
          std::get<Is>(models_).calc(local_data, x, u);
          data_collection.g_active.segment(g_offset<Is>(), DataType::G_SIZE) = local_data.g;
          data_collection.h_active.segment(h_offset<Is>(), DataType::H_SIZE) = local_data.h; }(), 0)...};
            }

            template <typename StateVector, typename ControlVector>
            void calcAllFastPath(const StateVector &x,
                                 const ControlVector &u,
                                 ConstraintDataCollection &data_collection)
            {
                calcAllFastPathImpl(x, u, data_collection, std::make_index_sequence<sizeof...(ConstraintModels)>{});
            }

            //-----------------------------------------------------------------
            // Cached path implementation:
            // Iterates over all tuple entries, and if the corresponding bit in active_mask_
            // is set, calls the model’s calc() and assigns its buffers into the output.
            template <std::size_t I = 0, typename StateVector, typename ControlVector>
            void calcAllCachedPathImpl(int &g_offset,
                                       int &h_offset,
                                       const StateVector &x,
                                       const ControlVector &u,
                                       ConstraintDataCollection &data_collection)
            {
                if constexpr (I < sizeof...(ConstraintModels))
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        using ModelType = std::tuple_element_t<I, std::tuple<ConstraintModels...>>;
                        using DataType = typename ModelType::DataType;
                        DataType local_data;
                        std::get<I>(models_).calc(local_data, x, u);
                        data_collection.g_active.segment(g_offset, DataType::G_SIZE) = local_data.g;
                        data_collection.h_active.segment(h_offset, DataType::H_SIZE) = local_data.h;
                        g_offset += DataType::G_SIZE;
                        h_offset += DataType::H_SIZE;
                    }
                    calcAllCachedPathImpl<I + 1>(g_offset, h_offset, x, u, data_collection);
                }
            }

            template <typename StateVector, typename ControlVector>
            void calcAllCachedPath(const StateVector &x,
                                   const ControlVector &u,
                                   ConstraintDataCollection &data_collection)
            {
                int g_offset = 0;
                int h_offset = 0;
                calcAllCachedPathImpl(g_offset, h_offset, x, u, data_collection);
            }

            //-----------------------------------------------------------------
            // Recursively sum the sizes of g and h for all active constraints.
            template <std::size_t I = 0>
            void updateDataCollectionSizesImpl(int &total_g, int &total_h) const
            {
                if constexpr (I < sizeof...(ConstraintModels))
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        using ModelType = std::tuple_element_t<I, std::tuple<ConstraintModels...>>;
                        using DataType = typename ModelType::DataType;
                        total_g += DataType::G_SIZE;
                        total_h += DataType::H_SIZE;
                    }
                    updateDataCollectionSizesImpl<I + 1>(total_g, total_h);
                }
            }

        }; // class ConstraintModelManager

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_manager_hpp__