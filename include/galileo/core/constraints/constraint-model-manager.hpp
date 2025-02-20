#ifndef __galileo_core_constraints_constraint_model_manager_hpp__
#define __galileo_core_constraints_constraint_model_manager_hpp__

#include "galileo/core/fwd.hpp"

#include <tuple>
#include <array>
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
        // - Before calling calc(), the user should call updateDataCollectionSizes()
        //   to resize h_active and g_active to exactly fit the sum of active constraints’ sizes.
        // - In calc(), if all constraints are active, a fast, fully unrolled code path is used;
        //   otherwise, a cached code path that tests each model’s active flag is used.
        template <typename... ConstraintModels_>
        class ConstraintModelManager
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ConstraintModels = std::tuple<ConstraintModels_...>;
            using ConstraintDatas = std::tuple<typename ConstraintModels_::ConstraintDataDerived...>;

            // full_mask_ is all ones (for the number of models).
            static constexpr std::size_t full_mask_ = (sizeof(ConstraintModels) >= 64 ? ~0ULL : ((1ULL << sizeof(ConstraintModels)) - 1));

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
                constexpr std::size_t idx = tuple_index<ConstraintType, ConstraintModels>::value;
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
                constexpr std::size_t idx = tuple_index<ConstraintType, ConstraintModels>::value;
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
            // buffers for all active constraints. Then resizes data accordingly.
            void updateDataCollectionSizes(ConstraintDataManager &data) const
            {
                int total_h = 0;
                int total_g = 0;
                updateDataCollectionSizesImpl<0>(total_h, total_g);
                data.h_active.resize(total_h);
                data.g_active.resize(total_g);

                data.h_eq_active.resize(total_h);
                data.g_lb_active.resize(total_g);
                data.g_ub_active.resize(total_g);

                data.Hx_active.resize(total_h, ConstraintModelManager::State_t::NDX);
                data.Hu_active.resize(total_h, ConstraintModelManager::ActuationModel_t::NU);
                data.Gx_active.resize(total_g, ConstraintModelManager::State_t::NDX);
                data.Gu_active.resize(total_g, ConstraintModelManager::ActuationModel_t::NU);
            }

            // calc:
            // Calls each active constraint’s calc() and assigns its output buffers into the
            // corresponding segments of data.
            // Two code paths are provided:
            //  - Fast path: when all constraints are active.
            //  - Cached path: when only a subset is active.
            template <typename StateVector, typename ControlVector>
            void calc(ConstraintDataManager &data,
                      const Eigen::MatrixBase<StateVector> &x,
                      const Eigen::MatrixBase<ControlVector> &u)
            {
                if (active_mask_ == full_mask_)
                {
                    // Fast path: all constraints active.
                    calcFastPath(data, x.derived(), u.derived());
                }
                else
                {
                    // Cached path: only some constraints active.
                    calcCachedPath(data, x.derived(), u.derived());
                }
            }

            // calcDiff:
            // Calls each active constraint’s calcDiff() and assigns its output buffers into the
            // corresponding segments of data.
            // Two code paths are provided:
            //  - Fast path: when all constraints are active.
            //  - Cached path: when only a subset is active.
            template <typename StateVector, typename ControlVector>
            void calcDiff(ConstraintDataManager &data,
                          const Eigen::MatrixBase<StateVector> &x,
                          const Eigen::MatrixBase<ControlVector> &u)
            {
                if (active_mask_ == full_mask_)
                {
                    // Fast path: all constraints active.
                    calcDiffFastPath(data, x.derived(), u.derived());
                }
                else
                {
                    // Cached path: only some constraints active.
                    calcDiffFastPath(data, x.derived(), u.derived());
                }
            }

        protected:
            //-----------------------------------------------------------------
            // Helpers to compute compile-time cumulative offsets.
            // For each index I, these return the starting offset (in h or g)
            // for the I-th constraint.
            template <std::size_t I>
            static constexpr int h_offset()
            {
                if constexpr (I == 0)
                    return 0;
                else
                    return h_offset<I - 1>() + std::tuple_element_t<I - 1, ConstraintModels>::ConstraintDataDerived::NH;
            }

            template <std::size_t I>
            static constexpr int g_offset()
            {
                if constexpr (I == 0)
                    return 0;
                else
                    return g_offset<I - 1>() + std::tuple_element_t<I - 1, ConstraintModels>::ConstraintDataDerived::NG;
            }

            //-----------------------------------------------------------------
            // Fast path implementation (all constraints active):
            // Unrolls the calls for each constraint using an index sequence.
            template <typename StateVector, typename ControlVector, std::size_t... Is>
            void calcFastPathImpl(ConstraintDataManager &data,
                                  const Eigen::MatrixBase<StateVector> &x,
                                  const Eigen::MatrixBase<ControlVector> &u,
                                  std::index_sequence<Is...>)
            {
                (void)std::initializer_list<int>{
                    ([&]()
                     {
                    using ConstraintModelDerived = std::tuple_element_t<Is, ConstraintModels>;
                    using ConstraintDataDerived = typename ConstraintModelDerived::ConstraintDataDerived;
                    ConstraintModelDerived local_model = std::get<Is>(models_);
                    ConstraintDataDerived local_data = std::get<ConstraintDataDerived>(data);

                    local_model.calc(local_data, x.derived(), u.derived());

                    data.h_active.template segment<ConstraintModelDerived::NH>(h_offset<Is>()) = local_data.h;
                    data.g_active.template segment<ConstraintModelDerived::NG>(g_offset<Is>()) = local_data.g;

                    data.h_eq_active.template segment<ConstraintModelDerived::NH>(h_offset<Is>()) = local_data.h_eq;
                    data.g_lb_active.template segment<ConstraintModelDerived::NG>(g_offset<Is>()) = local_data.g_lb;
                    data.g_ub_active.template segment<ConstraintModelDerived::NG>(g_offset<Is>()) = local_data.g_ub; }(), 0)...};
            }

            template <typename StateVector, typename ControlVector>
            void calcFastPath(ConstraintDataManager &data,
                              const Eigen::MatrixBase<StateVector> &x,
                              const Eigen::MatrixBase<ControlVector> &u)
            {
                calcFastPathImpl(data, x.derived(), u.derived(), std::make_index_sequence<sizeof(ConstraintModels)>{});
            }

            //-----------------------------------------------------------------
            // Cached path implementation:
            // Iterates over all tuple entries, and if the corresponding bit in active_mask_
            // is set, calls the model’s calc() and assigns its buffers into the output.
            template <std::size_t I = 0, typename StateVector, typename ControlVector>
            void calcCachedPathImpl(ConstraintDataManager &data,
                                    const Eigen::MatrixBase<StateVector> &x,
                                    const Eigen::MatrixBase<ControlVector> &u,
                                    int &h_offset,
                                    int &g_offset)
            {
                if constexpr (I < sizeof(ConstraintModels))
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        using ConstraintModelDerived = std::tuple_element_t<I, ConstraintModels>;
                        using ConstraintDataDerived = typename ConstraintModelDerived::ConstraintDataDerived;
                        ConstraintModelDerived local_model = std::get<I>(models_);
                        ConstraintDataDerived local_data = data.template get<ConstraintDataDerived>();

                        local_model.calc(local_data, x.derived(), u.derived());

                        data.h_active.template segment<ConstraintModelDerived::NH>(h_offset) = local_data.h;
                        data.g_active.template segment<ConstraintModelDerived::NG>(g_offset) = local_data.g;

                        data.h_eq_active.template segment<ConstraintModelDerived::NH>(h_offset) = local_data.h_eq;
                        data.g_lb_active.template segment<ConstraintModelDerived::NG>(g_offset) = local_data.g_lb;
                        data.g_ub_active.template segment<ConstraintModelDerived::NG>(g_offset) = local_data.g_ub;

                        h_offset += ConstraintModelDerived::NH;
                        g_offset += ConstraintModelDerived::NG;
                    }
                    calcCachedPathImpl<I + 1>(data, x.derived(), u.derived(), h_offset, g_offset);
                }
            }

            template <typename StateVector, typename ControlVector>
            void calcCachedPath(ConstraintDataManager &data,
                                const Eigen::MatrixBase<StateVector> &x,
                                const Eigen::MatrixBase<ControlVector> &u)
            {
                int h_offset = 0;
                int g_offset = 0;
                calcCachedPathImpl(data, x.derived(), u.derived(), h_offset, g_offset);
            }

            //-----------------------------------------------------------------
            // Recursively sum the sizes of g and h for all active constraints.
            template <std::size_t I = 0>
            void updateDataCollectionSizesImpl(int &total_h, int &total_g) const
            {
                if constexpr (I < sizeof(ConstraintModels))
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        using ConstraintModelDerived = std::tuple_element_t<I, ConstraintModels>;
                        using ConstraintDataDerived = typename ConstraintModelDerived::ConstraintDataDerived;
                        total_h += ConstraintModelDerived::NH;
                        total_g += ConstraintModelDerived::NG;
                    }
                    updateDataCollectionSizesImpl<I + 1>(total_h, total_g);
                }
            }

            //-----------------------------------------------------------------
            // Fast path implementation (all constraints active):
            // Unrolls the calls for each constraint using an index sequence.
            template <typename StateVector, typename ControlVector, std::size_t... Is>
            void calcDiffFastPathImpl(ConstraintDataManager &data,
                                      const Eigen::MatrixBase<StateVector> &x,
                                      const Eigen::MatrixBase<ControlVector> &u,
                                      std::index_sequence<Is...>)
            {
                (void)std::initializer_list<int>{
                    ([&]()
                     {
                    using ConstraintModelDerived = std::tuple_element_t<Is, ConstraintModels>;
                    using ConstraintDataDerived = typename ConstraintModelDerived::ConstraintDataDerived;
                    ConstraintModelDerived local_model = std::get<Is>(models_);
                    ConstraintDataDerived local_data = data.template get<ConstraintDataDerived>();

                    local_model.calcDiff(local_data, x.derived(), u.derived());

                    data.Hx_active.template block<ConstraintModelDerived::NH, ConstraintModelDerived::State_t::NDX>(h_offset<Is>(), 0) = local_data.Hx;
                    data.Hu_active.template block<ConstraintModelDerived::NH, ConstraintModelDerived::ActuationModel_t::NU>(h_offset<Is>(), 0) = local_data.Hu;
                    data.Gx_active.template block<ConstraintModelDerived::NG, ConstraintModelDerived::State_t::NDX>(g_offset<Is>(), 0) = local_data.Gx;
                    data.Gu_active.template block<ConstraintModelDerived::NG, ConstraintModelDerived::ActuationModel_t::NU>(g_offset<Is>(), 0) = local_data.Gu; }(), 0)...};
            }

            template <typename StateVector, typename ControlVector>
            void calcDiffFastPath(ConstraintDataManager &data,
                                  const Eigen::MatrixBase<StateVector> &x,
                                  const Eigen::MatrixBase<ControlVector> &u)
            {
                calcDiffFastPathImpl(data, x.derived(), u.derived(), std::make_index_sequence<sizeof(ConstraintModels)>{});
            }

            //-----------------------------------------------------------------
            // Cached path implementation:
            // Iterates over all tuple entries, and if the corresponding bit in active_mask_
            // is set, calls the model’s calc() and assigns its buffers into the output.
            template <std::size_t I = 0, typename StateVector, typename ControlVector>
            void calcDiffCachedPathImpl(ConstraintDataManager &data,
                                        const Eigen::MatrixBase<StateVector> &x,
                                        const Eigen::MatrixBase<ControlVector> &u,
                                        int &h_offset,
                                        int &g_offset)
            {
                if constexpr (I < sizeof(ConstraintModels))
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        using ConstraintModelDerived = std::tuple_element_t<I, ConstraintModels>;
                        using ConstraintDataDerived = typename ConstraintModelDerived::ConstraintDataDerived;
                        ConstraintModelDerived local_model = std::get<I>(models_);
                        ConstraintDataDerived local_data = data.template get<ConstraintDataDerived>();

                        local_model.calcDiff(local_data, x.derived(), u.derived());

                        data.Hx_active.template block<ConstraintModelDerived::NH, ConstraintModelDerived::State_t::NDX>(h_offset, 0) = local_data.Hx;
                        data.Hu_active.template block<ConstraintModelDerived::NH, ConstraintModelDerived::ActuationModel_t::NU>(h_offset, 0) = local_data.Hu;
                        data.Gx_active.template block<ConstraintModelDerived::NG, ConstraintModelDerived::State_t::NDX>(g_offset, 0) = local_data.Gx;
                        data.Gu_active.template block<ConstraintModelDerived::NG, ConstraintModelDerived::ActuationModel_t::NU>(g_offset, 0) = local_data.Gu;

                        h_offset += ConstraintModelDerived::NH;
                        g_offset += ConstraintModelDerived::NG;
                    }
                    calcDiffCachedPathImpl<I + 1>(data, x.derived(), u.derived(), h_offset, g_offset);
                }
            }

            template <typename StateVector, typename ControlVector>
            void calcDiffCachedPath(ConstraintDataManager &data,
                                    const Eigen::MatrixBase<StateVector> &x,
                                    const Eigen::MatrixBase<ControlVector> &u)
            {
                int h_offset = 0;
                int g_offset = 0;
                calcDiffCachedPathImpl(data, x.derived(), u.derived(), h_offset, g_offset);
            }

            // Tuple of constraint models.
            ConstraintModels models_;

            // Bit mask: bit i corresponds to the i-th constraint in the tuple.
            // A bit set to 1 means that constraint is active.
            std::size_t active_mask_;

        }; // class ConstraintModelManager

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_manager_hpp__