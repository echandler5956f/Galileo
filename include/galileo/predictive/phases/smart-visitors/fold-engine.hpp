#ifndef __galileo_predictive_phases_smart_visitors_fold_engine_hpp__
#define __galileo_predictive_phases_smart_visitors_fold_engine_hpp__

#include "galileo/predictive/phases/phase-generic.hpp"
#include "galileo/predictive/phases/smart-visitors/boundary-visitor.hpp"
#include "galileo/predictive/phases/smart-visitors/interior-visitor.hpp"
#include <vector>

namespace galileo
{
    namespace fusion
    {

        template <typename InteriorPropagator, typename BoundaryPropagator>
        concept IsFoldStateValid =
            std::is_same_v<typename InteriorPropagator::FoldStateType, typename BoundaryPropagator::FoldStateType>;

        template <typename InteriorPropagator, typename BoundaryPropagator>
        concept IsFoldReturnTypeValid =
            std::is_same_v<typename InteriorPropagator::ReturnType, typename BoundaryPropagator::ReturnType>;

        // Main sequential propagation engine with directional folding support
        template <typename InteriorPropagator,
                  typename BoundaryPropagator,
                  bool IsLeftFold_ = true> // true = forward propagation, false = backward propagation
            requires IsFoldStateValid<InteriorPropagator, BoundaryPropagator> &&
            IsFoldReturnTypeValid<InteriorPropagator, BoundaryPropagator>
        struct FoldEngineTpl
        {
        private:
            static constexpr bool IsLeftFold = IsLeftFold_;

            using InteriorBase = InteriorPropagatorBase<InteriorPropagator>;
            using BoundaryBase = BoundaryPropagatorBase<BoundaryPropagator>;

            using FoldStateType = typename traits<InteriorPropagator>::FoldStateType;
            using ReturnType = typename traits<InteriorPropagator>::ReturnType;

            // Directionally invariant means that regardless of the direction of the propagation,
            // the current phase is always considered the left phase, and the next (potentially reverse)
            // iterator is the right phase.
            static constexpr bool IsDirectionallyInvariant = traits<BoundaryPropagator>::IsDirectionallyInvariant;

        public:
            // Main propagation with args
            template <typename BasicSpec, template <typename> class CollectionTpl, typename ArgsTmp>
            static ReturnType run(const std::vector<PhaseModelTpl<BasicSpec, CollectionTpl>> &phase_models,
                                  std::vector<PhaseDataTpl<BasicSpec, CollectionTpl>> &phase_data,
                                  FoldStateType initial_state,
                                  ArgsTmp args)
            {
                return processPhases(phase_models, phase_data, initial_state, args);
            }

            // Main propagation without args
            template <typename BasicSpec, template <typename> class CollectionTpl>
            static ReturnType run(const std::vector<PhaseModelTpl<BasicSpec, CollectionTpl>> &phase_models,
                                  std::vector<PhaseDataTpl<BasicSpec, CollectionTpl>> &phase_data,
                                  FoldStateType initial_state)
            {
                return processPhases(phase_models, phase_data, initial_state);
            }

        private:
            // Unified phase processing with directional folding
            template <typename BasicSpec, template <typename> class CollectionTpl, typename ArgsTmp>
            static ReturnType processPhases(const std::vector<PhaseModelTpl<BasicSpec, CollectionTpl>> &phase_models,
                                            std::vector<PhaseDataTpl<BasicSpec, CollectionTpl>> &phase_data,
                                            FoldStateType initial_state,
                                            ArgsTmp args)
            {
                FoldStateType current_state = initial_state;

                const auto model_iterators = get_iterators(phase_models);
                auto data_iterators = get_iterators(phase_data);

                for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                     model_it != model_iterators.second;
                     ++model_it, ++data_it)
                {
                    current_state = processPhaseInterior(*model_it, *data_it, current_state, args);

                    // In forward iteration (IsLeftFold = true): next_model/data points to the subsequent phase
                    // In backward iteration (IsLeftFold = false): next_model/data points to the previous phase
                    const auto next_model = std::next(model_it);
                    auto next_data = std::next(data_it);
                    if (next_model != model_iterators.second)
                    {
                        if constexpr (IsDirectionallyInvariant || IsLeftFold)
                        {
                            current_state =
                                BoundaryBase::run(*model_it, *data_it, *next_model, *next_data, current_state, args);
                        }
                        else
                        {
                            current_state =
                                BoundaryBase::run(*next_model, *next_data, *model_it, *data_it, current_state, args);
                        }
                    }
                }

                return current_state;
            }

            template <typename BasicSpec, template <typename> class CollectionTpl>
            static ReturnType processPhases(const std::vector<PhaseModelTpl<BasicSpec, CollectionTpl>> &phase_models,
                                            std::vector<PhaseDataTpl<BasicSpec, CollectionTpl>> &phase_data,
                                            FoldStateType initial_state)
            {
                FoldStateType current_state = initial_state;

                const auto model_iterators = get_iterators(phase_models);
                auto data_iterators = get_iterators(phase_data);

                for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                     model_it != model_iterators.second;
                     ++model_it, ++data_it)
                {
                    current_state = processPhaseInterior(*model_it, *data_it, current_state);

                    // In forward iteration (IsLeftFold = true): next_model/data points to the subsequent phase
                    // In backward iteration (IsLeftFold = false): next_model/data points to the previous phase
                    const auto next_model = std::next(model_it);
                    auto next_data = std::next(data_it);
                    if (next_model != model_iterators.second)
                    {
                        if constexpr (IsDirectionallyInvariant || IsLeftFold)
                        {
                            current_state =
                                BoundaryBase::run(*model_it, *data_it, *next_model, *next_data, current_state);
                        }
                        else
                        {
                            current_state =
                                BoundaryBase::run(*next_model, *next_data, *model_it, *data_it, current_state);
                        }
                    }
                }

                return current_state;
            }

            // Process interior segments within a single phase
            template <typename BasicSpec, template <typename> class CollectionTpl, typename ArgsTmp>
            static ReturnType processPhaseInterior(const PhaseModelTpl<BasicSpec, CollectionTpl> &phase_model,
                                                   PhaseDataTpl<BasicSpec, CollectionTpl> &phase_data,
                                                   FoldStateType state,
                                                   ArgsTmp args)
            {
                InternalPhaseInteriorVisitor<BasicSpec, CollectionTpl, ArgsTmp> visitor(phase_data, state, args);
                return boost::apply_visitor(visitor, phase_model);
            }

            template <typename BasicSpec, template <typename> class CollectionTpl>
            static ReturnType processPhaseInterior(const PhaseModelTpl<BasicSpec, CollectionTpl> &phase_model,
                                                   PhaseDataTpl<BasicSpec, CollectionTpl> &phase_data,
                                                   FoldStateType state)
            {
                InternalPhaseInteriorVisitor<BasicSpec, CollectionTpl, NoArg> visitor(phase_data, state);
                return boost::apply_visitor(visitor, phase_model);
            }

            // Internal visitor for processing segments within a phase
            template <typename BasicSpec, template <typename> class CollectionTpl, typename ArgType>
            struct InternalPhaseInteriorVisitor : public boost::static_visitor<ReturnType>
            {
                using PhaseDataVariant_t = PhaseDataTpl<BasicSpec, CollectionTpl>;

                InternalPhaseInteriorVisitor(PhaseDataVariant_t &phase_data_, FoldStateType state_, ArgType args_)
                    : phase_data(phase_data_), state(state_), args(args_)
                {
                }

                template <typename PhaseModelType>
                ReturnType operator()(const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model) const
                {
                    using PhaseDataType = typename traits<PhaseModelType>::Data_t;

                    FoldStateType current_state = state;

                    const auto model_iterators = get_iterators(phase_model.derived().get_segments());
                    PhaseDataType &data = boost::get<PhaseDataType>(phase_data);
                    auto data_iterators = get_iterators(data.get_segments());

                    // In forward iteration (IsLeftFold = true): model_it/data_it points to the subsequent segment
                    // In backward iteration (IsLeftFold = false): model_it/data_it points to the previous segment
                    for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                         model_it != model_iterators.second;
                         ++model_it, ++data_it)
                    {
                        current_state = InteriorBase::run(*model_it, *data_it, current_state, args);
                    }

                    return current_state;
                }

                PhaseDataVariant_t &phase_data;
                FoldStateType state;
                ArgType args;
            };

            // Specialization for NoArg
            template <typename BasicSpec, template <typename> class CollectionTpl>
            struct InternalPhaseInteriorVisitor<BasicSpec, CollectionTpl, NoArg>
                : public boost::static_visitor<ReturnType>
            {
                using PhaseDataVariant_t = PhaseDataTpl<BasicSpec, CollectionTpl>;

                InternalPhaseInteriorVisitor(PhaseDataVariant_t &phase_data_, FoldStateType state_)
                    : phase_data(phase_data_), state(state_)
                {
                }

                template <typename PhaseModelType>
                ReturnType operator()(const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model) const
                {
                    using PhaseDataType = typename traits<PhaseModelType>::Data_t;

                    FoldStateType current_state = state;

                    const auto model_iterators = get_iterators(phase_model.derived().get_segments());
                    PhaseDataType &data = boost::get<PhaseDataType>(phase_data);
                    auto data_iterators = get_iterators(data.get_segments());

                    // In forward iteration (IsLeftFold = true): model_it/data_it points to the subsequent segment
                    // In backward iteration (IsLeftFold = false): model_it/data_it points to the previous segment
                    for (auto model_it = model_iterators.first, data_it = data_iterators.first;
                         model_it != model_iterators.second;
                         ++model_it, ++data_it)
                    {
                        current_state = InteriorBase::run(*model_it, *data_it, current_state);
                    }

                    return current_state;
                }

                PhaseDataVariant_t &phase_data;
                FoldStateType state;
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
        template <typename InteriorPropagator, typename BoundaryPropagator, bool IsLeftFold = true>
        using FoldTpl = FoldEngineTpl<InteriorPropagator, BoundaryPropagator, IsLeftFold>;

    } // namespace fusion

} // namespace galileo

#endif // __galileo_common_visitors_smart_visitors_fold_engine_hpp__
