#ifndef __galileo_predictive_segments_segment_data_erk_base_hpp__
#define __galileo_predictive_segments_segment_data_erk_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/segments/segment-model-erk-base.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        struct SegmentDataERKBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

        protected:
            inline SegmentDataERKBase()
            {
            }

            inline SegmentDataERKBase(const SegmentDataERKBase &clone)
            {
                *this = clone;
            }

            inline SegmentDataERKBase &operator=(const SegmentDataERKBase &clone)
            {
                return *this;
            }

        }; // struct SegmentDataERKBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_data_erk_base_hpp__
