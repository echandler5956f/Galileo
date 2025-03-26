#ifndef __galileo_core_costs_cost_data_manager_hpp__
#define __galileo_core_costs_cost_data_manager_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/core/costs/fwd.hpp"

namespace galileo
{
    namespace core
    {

        //-----------------------------------------------------------------
        // CostDataManager:
        // - Maintains aggregated vectors/matrices for *active* costs.
        // - Stores per-cost data in a tuple.
        // - Provides type-safe access to per-cost data.
        // - Used in conjunction with CostModelManager.
        template <typename... CostDatas_>
        class CostDataManager
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using CostDatas = std::tuple<CostDatas_...>;

            using VarScalar = typename std::tuple_element_t<0, CostDatas>::VarScalar;
            using NumScalar = typename std::tuple_element_t<0, CostDatas>::NumScalar;
            using Options = typename std::tuple_element_t<0, CostDatas>::Options;

            // Constructors
            template <typename... Datas>
            CostDataManager(Datas &&...datas)
                : data_(std::forward<Datas>(datas)...)
            {
            }

            // Accessors to the aggregated matrices
            auto &L() { return active_cost_matrices_.L; }
            auto &Lx() { return active_cost_matrices_.Lx; }
            auto &Lu() { return active_cost_matrices_.Lu; }
            auto &Lxx() { return active_cost_matrices_.Lxx; }
            auto &Lxu() { return active_cost_matrices_.Lxu; }
            auto &Luu() { return active_cost_matrices_.Luu; }

            // Public interface to retrieve typed data stored in the tuple.
            // This allows each cost model to access its specialized "CostDataDerived".
            template <typename T>
            T &get()
            {
                static_assert(HasType<T, DataTuple>::value,
                              "Requested type T not found in CostDataManager’s tuple.");
                return std::get<T>(data_);
            }

            template <typename T>
            const T &get() const
            {
                static_assert(HasType<T, DataTuple>::value,
                              "Requested type T not found in CostDataManager’s tuple.");
                return std::get<T>(data_);
            }

        private:
            // A simple helper metafunction to ensure T is in the tuple.
            template <typename T, typename Tuple>
            struct HasType;

            template <typename T>
            struct HasType<T, std::tuple<>> : std::false_type
            {
            };

            template <typename T, typename... Ts>
            struct HasType<T, std::tuple<T, Ts...>> : std::true_type
            {
            };

            template <typename T, typename U, typename... Ts>
            struct HasType<T, std::tuple<U, Ts...>> : HasType<T, std::tuple<Ts...>>
            {
            };

            // Aggregated outputs for the "active" costs
            CostMatricesDynamic<VarScalar, NumScalar, Options> active_cost_matrices_;

            // Storage of per-cost data
            CostDatas data_;

        }; // class CostDataManager

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_data_manager_hpp__
