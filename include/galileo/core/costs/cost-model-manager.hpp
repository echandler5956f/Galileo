#ifndef __galileo_core_costs_cost_model_manager_hpp__
#define __galileo_core_costs_cost_model_manager_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/core/costs/cost-data-manager.hpp"

#include <tuple>
#include <array>
#include <initializer_list>

namespace galileo
{

    namespace core
    {
        //-----------------------------------------------------------------
        // CostModelManager:
        // - Holds a tuple of cost models and an active bitmask.
        // - Provides type–safe activate()/deactivate() functions.
        // - In calc(), a cached code path that tests each model’s active flag is used.
        template <typename... CostModels_>
        class CostModelManager
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using CostModels = std::tuple<CostModels_...>;
            using CostDatas = std::tuple<typename CostModels_::CostDataDerived...>;

            using NumScalar = typename std::tuple_element_t<0, CostModels>::NumScalar;

            static constexpr std::size_t N = std::tuple_size<CostModels>::value;

            // full_mask_ is all ones (for the number of models).
            static constexpr std::size_t full_mask_ = (N >= 64 ? ~0ULL : ((1ULL << N) - 1));

            // Constructor: all costs active by default.
            template <typename... Models>
            CostModelManager(Models &&...models)
                : models_(std::forward<Models>(models)...),
                  active_mask_(full_mask_)
            {
            }

            // --- Single-cost activation:
            template <typename CostType>
            void activate()
            {
                constexpr std::size_t idx = tuple_index<CostType, CostModels>::value;
                active_mask_ |= (1ULL << idx);
            }

            // Variadic activation (only enabled when at least 2 types are provided):
            template <typename First, typename Second, typename... Others>
            void activate()
            {
                activate<First>();
                activate<Second, Others...>();
            }

            // --- Single-cost deactivation:
            template <typename CostType>
            void deactivate()
            {
                constexpr std::size_t idx = tuple_index<CostType, CostModels>::value;
                active_mask_ &= ~(1ULL << idx);
            }

            // Variadic deactivation (only enabled when at least 2 types are provided):
            template <typename First, typename Second, typename... Others>
            void deactivate()
            {
                deactivate<First>();
                deactivate<Second, Others...>();
            }

            // calc:
            // Calls each active cost’s calc() and assigns its output buffers into the
            // corresponding segments of data.
            // Cached path: when only a subset is active.
            template <typename StateVector, typename ControlVector>
            void calc(CostDataManager &data,
                      const Eigen::MatrixBase<StateVector> &x,
                      const Eigen::MatrixBase<ControlVector> &u)
            {
                data.L() = NumScalar(0.);

                // Cached path: only some costs are active.
                calcCachedPathImpl(data, x.derived(), u.derived());
            }

            // calcDiff:
            // Calls each active cost’s calcDiff() and assigns its output buffers into the
            // corresponding segments of data.
            // Cached path: when only a subset is active.
            template <typename StateVector, typename ControlVector>
            void calcDiff(CostDataManager &data,
                          const Eigen::MatrixBase<StateVector> &x,
                          const Eigen::MatrixBase<ControlVector> &u)
            {
                data.Lx().setZero();
                data.Lu().setZero();
                data.Lxx().setZero();
                data.Lxu().setZero();
                data.Luu().setZero();

                // Cached path: only some costs are active.
                calcDiffCachedPath(data, x.derived(), u.derived());
            }

        protected:
            //-----------------------------------------------------------------
            // Cached path implementation:
            // Iterates over all tuple entries, and if the corresponding bit in active_mask_
            // is set, calls the model’s calc() and assigns its buffers into the output.
            template <std::size_t I = 0, typename StateVector, typename ControlVector>
            void calcCachedPath(CostDataManager &data,
                                const Eigen::MatrixBase<StateVector> &x,
                                const Eigen::MatrixBase<ControlVector> &u)
            {
                if constexpr (I < N)
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        auto &m_i = std::get<I>(models_);
                        auto &d_i = data.template get<typename decltype(m_i)::CostDataDerived>();

                        m_i.calc(d_i, x.derived(), u.derived());

                        data.L() += m_i.weight() * d_i.L();
                    }
                    calcCachedPath<I + 1>(data, x.derived(), u.derived());
                }
            }

            //-----------------------------------------------------------------
            // Cached path implementation:
            // Iterates over all tuple entries, and if the corresponding bit in active_mask_
            // is set, calls the model’s calc() and assigns its buffers into the output.
            template <std::size_t I = 0, typename StateVector, typename ControlVector>
            void calcDiffCachedPath(CostDataManager &data,
                                    const Eigen::MatrixBase<StateVector> &x,
                                    const Eigen::MatrixBase<ControlVector> &u)
            {
                if constexpr (I < N)
                {
                    if (active_mask_ & (1ULL << I))
                    {
                        auto &m_i = std::get<I>(models_);
                        auto &d_i = data.template get<typename decltype(m_i)::CostDataDerived>();

                        m_i.calcDiff(d_i, x.derived(), u.derived());

                        data.Lx() += m_i.weight() * d_i.Lx();
                        data.Lu() += m_i.weight() * d_i.Lu();
                        data.Lxx() += m_i.weight() * d_i.Lxx();
                        data.Lxu() += m_i.weight() * d_i.Lxu();
                        data.Luu() += m_i.weight() * d_i.Luu();
                    }
                    calcDiffCachedPath<I + 1>(data, x.derived(), u.derived());
                }
            }

            // Tuple of cost models.
            CostModels models_;

            // Bit mask: bit i corresponds to the i-th cost in the tuple.
            // A bit set to 1 means that cost is active.
            std::size_t active_mask_;

        }; // class CostModelManager

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_cost_model_manager_hpp__