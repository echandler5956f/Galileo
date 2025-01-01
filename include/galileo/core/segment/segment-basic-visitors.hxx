#ifndef __galileo_core_segment_basic_visitors_hxx__
#define __galileo_core_segment_basic_visitors_hxx__

#include "galileo/core/segment/segment-basic-visitors.hpp"
#include "galileo/core/visitor.hpp"

namespace galileo
{

  template <typename Scalar, int Options, template <typename S, int O> class SegmentCollectionTpl>
  struct CreateSegmentData : boost::static_visitor<SegmentDataTpl<Scalar, Options, SegmentCollectionTpl>>
  {
    typedef SegmentCollectionTpl<Scalar, Options> SegmentCollection;
    typedef typename SegmentCollection::SegmentModelVariant SegmentModelVariant;
    typedef SegmentDataTpl<Scalar, Options, SegmentCollectionTpl> SegmentDataVariant;

    template <typename SegmentModelDerived>
    SegmentDataVariant operator()(const SegmentModelBase<SegmentModelDerived> &segment_model) const
    {
      return SegmentDataVariant(segment_model.createData());
    }

    static SegmentDataVariant run(const SegmentModelVariant &segment_model)
    {
      return boost::apply_visitor(CreateSegmentData(), segment_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class SegmentCollectionTpl>
  inline SegmentDataTpl<Scalar, Options, SegmentCollectionTpl>
  createData(const SegmentModelTpl<Scalar, Options, SegmentCollectionTpl> &segment_model)
  {
    return CreateSegmentData<Scalar, Options, SegmentCollectionTpl>::run(segment_model);
  }

} // namespace galileo

#endif // __galileo_core_segment_basic_visitors_hxx__
