#ifndef __galileo_predictive_phases_fold_visitors_fold_engine_hpp__
#define __galileo_predictive_phases_fold_visitors_fold_engine_hpp__

#include "galileo/predictive/phases/fold-visitors/boundary-visitor.hpp"
#include "galileo/predictive/phases/fold-visitors/interior-visitor.hpp"
#include "galileo/predictive/phases/phase-generic.hpp"
#include <vector>

namespace galileo
{
    namespace fusion
    {

        // Main sequential propagation engine with directional folding support
        template <
            typename StateType,
            typename InteriorPropagator,
            typename BoundaryPropagator,
            bool IsLeftFold = true, // true = forward propagation, false = backward propagation
            typename ReturnType = StateType>
        struct FoldEngineTpl
        {
        private:
            using InteriorBase = InteriorPropagatorBase<StateType, InteriorPropagator, ReturnType>;
            using BoundaryBase = BoundaryPropagatorBase<StateType, BoundaryPropagator, ReturnType>;

        public:
            // Main propagation with args
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const std::vector<PhaseModelTpl<PhaseSpec, CollectionTpl>> &phase_models,
                std::vector<PhaseDataTpl<PhaseSpec, CollectionTpl>> &phase_data,
                StateType initial_state,
                ArgsTmp args)
            {
                return processPhases(phase_models, phase_data, initial_state, args);
            }

            // Main propagation without args
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl>
            static ReturnType run(
                const std::vector<PhaseModelTpl<PhaseSpec, CollectionTpl>> &phase_models,
                std::vector<PhaseDataTpl<PhaseSpec, CollectionTpl>> &phase_data,
                StateType initial_state)
            {
                return processPhases(phase_models, phase_data, initial_state);
            }

        private:
            // Unified phase processing with directional folding
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType processPhases(
                const std::vector<PhaseModelTpl<PhaseSpec, CollectionTpl>> &phase_models,
                std::vector<PhaseDataTpl<PhaseSpec, CollectionTpl>> &phase_data,
                StateType initial_state,
                ArgsTmp args)
            {
                StateType current_state = initial_state;

                auto model_iterators = get_iterators(phase_models);
                auto data_iterators = get_iterators(phase_data);

                for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                     model_it != model_iterators.second; ++model_it, ++data_it)
                {
                    current_state = processPhaseInterior(*model_it, *data_it, current_state, args);

                    // In forward iteration (IsLeftFold = true): next_model/data points to the subsequent phase
                    // In backward iteration (IsLeftFold = false): next_model/data points to the previous phase
                    auto next_model = std::next(model_it);
                    auto next_data = std::next(data_it);
                    if (next_model != model_iterators.second)
                    {
                        current_state = BoundaryBase::run(*model_it, *data_it, *next_model, *next_data, current_state, args);
                    }
                }

                return current_state;
            }

            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl>
            static ReturnType processPhases(
                const std::vector<PhaseModelTpl<PhaseSpec, CollectionTpl>> &phase_models,
                std::vector<PhaseDataTpl<PhaseSpec, CollectionTpl>> &phase_data,
                StateType initial_state)
            {
                StateType current_state = initial_state;

                auto model_iterators = get_iterators(phase_models);
                auto data_iterators = get_iterators(phase_data);

                for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                     model_it != model_iterators.second; ++model_it, ++data_it)
                {
                    current_state = processPhaseInterior(*model_it, *data_it, current_state);

                    // In forward iteration (IsLeftFold = true): next_model/data points to the subsequent phase
                    // In backward iteration (IsLeftFold = false): next_model/data points to the previous phase
                    auto next_model = std::next(model_it);
                    auto next_data = std::next(data_it);
                    if (next_model != model_iterators.second)
                    {
                        current_state = BoundaryBase::run(*model_it, *data_it, *next_model, *next_data, current_state);
                    }
                }

                return current_state;
            }

            // Process interior segments within a single phase
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType processPhaseInterior(
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &phase_model,
                PhaseDataTpl<PhaseSpec, CollectionTpl> &phase_data,
                StateType state,
                ArgsTmp args)
            {
                InternalPhaseInteriorVisitor<PhaseSpec, CollectionTpl, ArgsTmp> visitor(phase_data, state, args);
                return boost::apply_visitor(visitor, phase_model);
            }

            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl>
            static ReturnType processPhaseInterior(
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &phase_model,
                PhaseDataTpl<PhaseSpec, CollectionTpl> &phase_data,
                StateType state)
            {
                InternalPhaseInteriorVisitor<PhaseSpec, CollectionTpl, NoArg> visitor(phase_data, state);
                return boost::apply_visitor(visitor, phase_model);
            }

            // Internal visitor for processing segments within a phase
            template <typename PhaseSpec, template <typename> class CollectionTpl, typename ArgType>
            struct InternalPhaseInteriorVisitor : public boost::static_visitor<ReturnType>
            {
                using PhaseDataVariant_t = PhaseDataTpl<PhaseSpec, CollectionTpl>;

                InternalPhaseInteriorVisitor(PhaseDataVariant_t &phase_data_, StateType state_, ArgType args_)
                    : phase_data(phase_data_), state(state_), args(args_) {}

                template <typename PhaseModelType>
                ReturnType operator()(const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model) const
                {
                    using PhaseDataType = typename traits<PhaseModelType>::Data_t;

                    StateType current_state = state;

                    auto model_iterators = get_iterators(phase_model.derived().get_segments_models());
                    PhaseDataType &data = boost::get<PhaseDataType>(phase_data);
                    auto data_iterators = get_iterators(data.get_segments_data());

                    // In forward iteration (IsLeftFold = true): model_it/data_it points to the subsequent segment
                    // In backward iteration (IsLeftFold = false): model_it/data_it points to the previous segment
                    for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                         model_it != model_iterators.second; ++model_it, ++data_it)
                    {
                        current_state = InteriorBase::run(*model_it, *data_it, current_state, args);
                    }

                    return current_state;
                }

                PhaseDataVariant_t &phase_data;
                StateType state;
                ArgType args;
            };

            // Specialization for NoArg
            template <typename PhaseSpec, template <typename> class CollectionTpl>
            struct InternalPhaseInteriorVisitor<PhaseSpec, CollectionTpl, NoArg> : public boost::static_visitor<ReturnType>
            {
                using PhaseDataVariant_t = PhaseDataTpl<PhaseSpec, CollectionTpl>;

                InternalPhaseInteriorVisitor(PhaseDataVariant_t &phase_data_, StateType state_)
                    : phase_data(phase_data_), state(state_) {}

                template <typename PhaseModelType>
                ReturnType operator()(const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model) const
                {
                    using PhaseDataType = typename traits<PhaseModelType>::Data_t;

                    StateType current_state = state;

                    auto model_iterators = get_iterators(phase_model.derived().get_segments_models());
                    PhaseDataType &data = boost::get<PhaseDataType>(phase_data);
                    auto data_iterators = get_iterators(data.get_segments_data());

                    // In forward iteration (IsLeftFold = true): model_it/data_it points to the subsequent segment
                    // In backward iteration (IsLeftFold = false): model_it/data_it points to the previous segment
                    for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                         model_it != model_iterators.second; ++model_it, ++data_it)
                    {
                        current_state = InteriorBase::run(*model_it, *data_it, current_state);
                    }

                    return current_state;
                }

                PhaseDataVariant_t &phase_data;
                StateType state;
            };

            // Get iterators for the container
            // Uses universal reference to avoid duplicating const vs non-const versions
            template <typename Container>
            static auto get_iterators(Container &&container)
            {
                if constexpr (IsLeftFold)
                {
                    return std::make_pair(container.begin(), container.end());
                }
                else
                {
                    return std::make_pair(container.rbegin(), container.rend());
                }
            }

        }; // struct FoldEngineTpl

        // Convenience alias matching the expected user interface
        template <
            typename StateType,
            typename InteriorPropagator,
            typename BoundaryPropagator,
            bool IsLeftFold = true>
        using FoldTpl = FoldEngineTpl<StateType, InteriorPropagator, BoundaryPropagator, IsLeftFold, StateType>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_common_visitors_fold_visitors_fold_engine_hpp__
