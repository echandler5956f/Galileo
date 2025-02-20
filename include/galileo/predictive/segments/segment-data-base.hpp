#ifndef __galileo_predictive_segments_segment_data_base_hpp__
#define __galileo_predictive_segments_segment_data_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/segments/segment-model-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        struct SegmentDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using NodeDerived = typename traits<Derived>::NodeDerived;
            using SegmentDerived = typename traits<Derived>::SegmentDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_CONSTANTS(SegmentDerived);

            GALILEO_NODE_DATA_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_DATA_TYPEDEF(SegmentDerived);

            NodeDataVector nodes;
            ControlParamData_t *control;
            C_t C;
            Ck_t Ck;
            Cw_t Cw;

        protected:
            inline SegmentDataBase()
            {
            }

            inline SegmentDataBase(const SegmentDataBase &clone)
            {
                *this = clone;
            }

            inline SegmentDataBase &operator=(const SegmentDataBase &clone)
            {
                return *this;
            }

        }; // struct SegmentDataBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_data_base_hpp__
