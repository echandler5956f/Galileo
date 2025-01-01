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
        typedef _Scalar Scalar;
        enum
        {
            Options = _Options
        };

        typedef boost::variant<
            SegmentModelVoid>
            SegmentModelVariant;

        typedef boost::variant<
            SegmentDataVoid>
            SegmentDataVariant;
    };

    typedef SegmentCollectionDefault::SegmentModelVariant SegmentModelVariant;
    typedef SegmentCollectionDefault::SegmentDataVariant SegmentDataVariant;

} // namespace galileo

#endif // __galileo_core_segment_collections_hpp__
