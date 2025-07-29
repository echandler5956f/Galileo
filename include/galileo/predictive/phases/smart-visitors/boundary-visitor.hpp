#ifndef __galileo_predictive_phases_smart_visitors_boundary_visitor_hpp__
#define __galileo_predictive_phases_smart_visitors_boundary_visitor_hpp__

#include "galileo/common/visitors/binary-visitor.hpp"
#include "galileo/predictive/phases/phase-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Boundary propagator base for phase-to-phase transitions (model-model pairs)
        template <typename VisitorDerived_>
        struct BoundaryPropagatorBase
        {
        private:
            template <typename LeftPhaseModelType>
            using LeftPhaseModelBaseOf_t = PhaseModelBase<LeftPhaseModelType, typename traits<LeftPhaseModelType>::PS>;
            template <typename RightPhaseModelType>
            using RightPhaseModelBaseOf_t = PhaseModelBase<RightPhaseModelType, typename traits<RightPhaseModelType>::PS>;

            template <typename LeftPhaseDataType>
            using LeftPhaseDataBaseOf_t = PhaseDataBase<LeftPhaseDataType, typename traits<LeftPhaseDataType>::PS>;
            template <typename RightPhaseDataType>
            using RightPhaseDataBaseOf_t = PhaseDataBase<RightPhaseDataType, typename traits<RightPhaseDataType>::PS>;

        public:
            using VisitorDerived = VisitorDerived_;
            using FoldStateType = typename traits<VisitorDerived>::FoldStateType;
            static constexpr bool IsDirectionallyInvariant = traits<VisitorDerived>::IsDirectionallyInvariant;
            using ReturnType = typename traits<VisitorDerived>::ReturnType;

            // Phase boundary transformation with args
            template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData, typename ArgsTmp>
            static ReturnType run(
                const LeftPhaseModelBaseOf_t<CurrentPhaseModel> &current_phase_model,
                const LeftPhaseDataBaseOf_t<CurrentPhaseData> &current_phase_data,
                const RightPhaseModelBaseOf_t<NextPhaseModel> &next_phase_model,
                const RightPhaseDataBaseOf_t<NextPhaseData> &next_phase_data,
                FoldStateType state,
                ArgsTmp args)
            {
                return bf::invoke(
                    &VisitorDerived::template algo<CurrentPhaseModel, CurrentPhaseData, NextPhaseModel, NextPhaseData>,
                    gf::append(
                        boost::ref(current_phase_model.derived()),
                        boost::ref(current_phase_data.derived()),
                        boost::ref(next_phase_model.derived()),
                        boost::ref(next_phase_data.derived()),
                        state,
                        args));
            }

            // Phase boundary transformation without args
            template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData>
            static ReturnType run(
                const LeftPhaseModelBaseOf_t<CurrentPhaseModel> &current_phase_model,
                const LeftPhaseDataBaseOf_t<CurrentPhaseData> &current_phase_data,
                const RightPhaseModelBaseOf_t<NextPhaseModel> &next_phase_model,
                const RightPhaseDataBaseOf_t<NextPhaseData> &next_phase_data,
                FoldStateType state)
            {
                return VisitorDerived::template algo<CurrentPhaseModel, CurrentPhaseData, NextPhaseModel, NextPhaseData>(
                    current_phase_model.derived(), current_phase_data.derived(), next_phase_model.derived(), next_phase_data.derived(), state);
            }

            // Binary visitor dispatch for type-erased phases
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl,
                typename ArgsTmp>
            static ReturnType run(
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &current_phase_model,
                const PhaseDataTpl<PhaseSpec, CollectionTpl> &current_phase_data,
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &next_phase_model,
                const PhaseDataTpl<PhaseSpec, CollectionTpl> &next_phase_data,
                FoldStateType state,
                ArgsTmp args)
            {
                InternalBoundaryVisitor<PhaseSpec, CollectionTpl, ArgsTmp> visitor(state, args);
                return boost::apply_visitor(visitor, current_phase_model, current_phase_data, next_phase_model, next_phase_data);
            }

            // Binary visitor dispatch without args
            template <
                typename PhaseSpec,
                template <typename> class CollectionTpl>
            static ReturnType run(
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &current_phase_model,
                const PhaseDataTpl<PhaseSpec, CollectionTpl> &current_phase_data,
                const PhaseModelTpl<PhaseSpec, CollectionTpl> &next_phase_model,
                const PhaseDataTpl<PhaseSpec, CollectionTpl> &next_phase_data,
                FoldStateType state)
            {
                InternalBoundaryVisitor<PhaseSpec, CollectionTpl, NoArg> visitor(state);
                return boost::apply_visitor(visitor, current_phase_model, current_phase_data, next_phase_model, next_phase_data);
            }

        private:
            // Internal visitor for binary dispatch on phase boundaries
            template <typename PhaseSpec, template <typename> class CollectionTpl, typename ArgType>
            struct InternalBoundaryVisitor : public boost::static_visitor<ReturnType>
            {
                InternalBoundaryVisitor(FoldStateType state_, ArgType args_)
                    : state(state_), args(args_) {}

                template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData>
                ReturnType operator()(const LeftPhaseModelBaseOf_t<CurrentPhaseModel> &current_phase_model,
                                      const LeftPhaseDataBaseOf_t<CurrentPhaseData> &current_phase_data,
                                      const RightPhaseModelBaseOf_t<NextPhaseModel> &next_phase_model,
                                      const RightPhaseDataBaseOf_t<NextPhaseData> &next_phase_data) const
                {
                    return BoundaryPropagatorBase::run(current_phase_model, current_phase_data, next_phase_model, next_phase_data, state, args);
                }

                FoldStateType state;
                ArgType args;
            };

            // Specialization for NoArg
            template <typename PhaseSpec, template <typename> class CollectionTpl>
            struct InternalBoundaryVisitor<PhaseSpec, CollectionTpl, NoArg> : public boost::static_visitor<ReturnType>
            {
                InternalBoundaryVisitor(FoldStateType state_) : state(state_) {}

                template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData>
                ReturnType operator()(const LeftPhaseModelBaseOf_t<CurrentPhaseModel> &current_phase_model,
                                      const LeftPhaseDataBaseOf_t<CurrentPhaseData> &current_phase_data,
                                      const RightPhaseModelBaseOf_t<NextPhaseModel> &next_phase_model,
                                      const RightPhaseDataBaseOf_t<NextPhaseData> &next_phase_data) const
                {
                    return BoundaryPropagatorBase::run(current_phase_model, current_phase_data, next_phase_model, next_phase_data, state);
                }

                FoldStateType state;
            };

        }; // struct BoundaryPropagatorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_predictive_phases_smart_visitors_boundary_visitor_hpp__
