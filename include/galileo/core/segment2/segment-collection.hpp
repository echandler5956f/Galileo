#ifndef __galileo_core_segment_collections_hpp__
#define __galileo_core_segment_collections_hpp__

#include "galileo/core/segment/fwd.hpp"
// #include "galileo/core/segment/segments.hpp"

#include <boost/variant.hpp>

namespace galileo
{

    template <typename _Scalar, int _Options>
    struct SegmentCollectionDefaultTpl
    {
        using Scalar = _Scalar;
        enum
        {
            Options = _Options
        };

        using SegmentModelVariant = boost::variant<SegmentModelVoid>;
        using SegmentDataVariant = boost::variant<SegmentDataVoid>;
    };

    using SegmentModelVariant = typename SegmentCollectionDefault::SegmentModelVariant;
    using SegmentDataVariant = typename SegmentCollectionDefault::SegmentDataVariant;

} // namespace galileo

#endif // __galileo_core_segment_collections_hpp__
