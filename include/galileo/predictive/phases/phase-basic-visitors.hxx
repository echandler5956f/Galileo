#ifndef __galileo_predictive_phases_phase_basic_visitors_hxx__
#define __galileo_predictive_phases_phase_basic_visitors_hxx__

#include "galileo/predictive/phases/phase-basic-visitors.hpp"
#include "galileo/predictive/visitor/phase-unary-visitor.hpp"

namespace galileo
{

    namespace predictive
    {

        struct CalcSegmentFromPhaseVisitor
            : fusion::PhaseUnaryVisitorBase<CalcSegmentFromPhaseVisitor>
        {
            typedef boost::fusion::vector<std::size_t> ArgsType;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const std::size_t &index)
            {
                phase_model.segments_[index].calc(phase_data[index]);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline void calc_diff_segment_from_phase(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &index)
        {
            typedef CalcSegmentFromPhaseVisitor Algo;

            Algo::run(phase_model, phase_data, typename Algo::ArgsType(index));
        }

        struct CalcDiffSegmentFromPhaseVisitor
            : fusion::PhaseUnaryVisitorBase<CalcDiffSegmentFromPhaseVisitor>
        {
            typedef boost::fusion::vector<std::size_t> ArgsType;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const std::size_t &index)
            {
                phase_model.segments_[index].calcDiff(phase_data[index]);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline void calc_diff_segment_from_phase(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &index)
        {
            typedef CalcDiffSegmentFromPhaseVisitor Algo;

            Algo::run(phase_model, phase_data, typename Algo::ArgsType(index));
        }

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hxx__