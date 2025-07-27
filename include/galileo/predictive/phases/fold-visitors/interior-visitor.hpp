#ifndef __galileo_predictive_phases_fold_visitors_interior_visitor_hpp__
#define __galileo_predictive_phases_fold_visitors_interior_visitor_hpp__

#include "galileo/common/visitors/unary-visitor.hpp"
#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/segments/segment-erk-base.hpp"

namespace galileo
{
    namespace fusion
    {

        // Interior propagator base for segment-level transformations (model-data pairs)
        template <typename FoldStateType, typename VisitorDerived, typename ReturnType = FoldStateType>
        struct InteriorPropagatorBase
        {
        private:
            template <typename SegmentModelType>
            using SegmentModelBaseOf_t = SegmentERKModelBase<SegmentModelType, typename traits<SegmentModelType>::PS>;

            template <typename SegmentDataType>
            using SegmentDataBaseOf_t = SegmentERKDataBase<SegmentDataType, typename traits<SegmentDataType>::PS>;

        public:
            // Segment model-data transformation with args
            template <typename SegmentModel, typename SegmentData, typename ArgsTmp>
            static ReturnType run(
                const SegmentModelBaseOf_t<SegmentModel> &segment_model,
                SegmentDataBaseOf_t<SegmentData> &segment_data,
                FoldStateType state,
                ArgsTmp args)
            {
                return bf::invoke(
                    &VisitorDerived::template algo<SegmentModel, SegmentData>,
                    gf::append(
                        boost::ref(segment_model.derived()),
                        boost::ref(segment_data.derived()),
                        state,
                        args));
            }

            // Segment model-data transformation without args
            template <typename SegmentModel, typename SegmentData>
            static ReturnType run(
                const SegmentModelBaseOf_t<SegmentModel> &segment_model,
                SegmentDataBaseOf_t<SegmentData> &segment_data,
                FoldStateType state)
            {
                return VisitorDerived::template algo<SegmentModel, SegmentData>(segment_model.derived(), segment_data.derived(), state);
            }

        }; // struct InteriorPropagatorBase

    } // namespace fusion

} // namespace galileo

#endif // __galileo_predictive_phases_fold_visitors_interior_visitor_hpp__
