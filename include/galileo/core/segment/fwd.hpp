#ifndef __galileo_core_segment_fwd_hpp__
#define __galileo_core_segment_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    struct SegmentModelVoid
    {
    };

    struct SegmentDataVoid
    {
    };

    template <typename Scalar, int Options = context::Options>
    struct SegmentCollectionDefaultTpl;
    using SegmentCollectionDefault = SegmentCollectionDefaultTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class SegmentCollectionTpl = SegmentCollectionDefaultTpl>
    struct SegmentModelTpl;
    using SegmentModel = SegmentModelTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class SegmentCollectionTpl = SegmentCollectionDefaultTpl>
    struct SegmentDataTpl;
    using SegmentData = SegmentDataTpl<context::Scalar>;

} // namespace galileo

#include "galileo/core/fwd.hpp"

#endif // __galileo_core_segment_fwd_hpp__
