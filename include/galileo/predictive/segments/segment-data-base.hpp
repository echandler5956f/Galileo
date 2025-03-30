#ifndef __galileo_predictive_segments_segment_data_base_hpp__
#define __galileo_predictive_segments_segment_data_base_hpp__

#include "galileo/predictive/segments/segment-base.hpp"
#include "galileo/predictive/segments/segment-model-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename Derived, typename PhaseSpec>
        struct SegmentDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            // NodeDataVector nodes;
            // ControlParamData_t *control;
            // S_t S;
            // Sk_t Sk;
            // Sw_t Sw;

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
