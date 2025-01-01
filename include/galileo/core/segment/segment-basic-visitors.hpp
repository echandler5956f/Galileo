#ifndef __galileo_core_segment_basic_visitors_hpp__
#define __galileo_core_segment_basic_visitors_hpp__

#include "galileo/core/segment/fwd.hpp"

namespace galileo
{

    template <typename Scalar, int Options, template <typename S, int O> class SegmentCollectionTpl>
    inline SegmentDataTpl<Scalar, Options, SegmentCollectionTpl>
    createData(const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model);

} // namespace galileo

/* --- Details -------------------------------------------------------------------- */
// Included later
// #include "galileo/core/segment/core-basic-visitors.hxx"

#endif // __galileo_core_segment_basic_visitors_hpp__
