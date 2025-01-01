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
    typedef SegmentCollectionDefaultTpl<context::Scalar> SegmentCollectionDefault;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class SegmentCollectionTpl = SegmentCollectionDefaultTpl>
    struct SegmentModelTpl;
    typedef SegmentModelTpl<context::Scalar> SegmentModel;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class SegmentCollectionTpl = SegmentCollectionDefaultTpl>
    struct SegmentDataTpl;
    typedef SegmentDataTpl<context::Scalar> SegmentData;

} // namespace galileo

#include "galileo/core/fwd.hpp"

#endif // __galileo_core_segment_fwd_hpp__
