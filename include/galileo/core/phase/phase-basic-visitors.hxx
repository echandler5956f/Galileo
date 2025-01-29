#ifndef __galileo_core_phase_basic_visitors_hxx__
#define __galileo_core_phase_basic_visitors_hxx__

#include "galileo/core/phase/phase-basic-visitors.hpp"
#include "galileo/core/visitor.hpp"

namespace galileo
{

  template <typename Scalar, int Options, template <typename S, int O> class PhaseCollectionTpl>
  struct CreatePhaseData : boost::static_visitor<PhaseDataTpl<Scalar, Options, PhaseCollectionTpl>>
  {
    using PhaseCollection = PhaseCollectionTpl<Scalar, Options>;
    using PhaseModelVariant = typename PhaseCollection::PhaseModelVariant;
    using PhaseDataVariant = PhaseDataTpl<Scalar, Options, PhaseCollectionTpl>;

    template <typename PhaseModelDerived>
    PhaseDataVariant operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
    {
      return PhaseDataVariant(phase_model.createData());
    }

    static PhaseDataVariant run(const PhaseModelVariant &phase_model)
    {
      return boost::apply_visitor(CreatePhaseData(), phase_model);
    }
  };

  template <typename Scalar, int Options, template <typename S, int O> class PhaseCollectionTpl>
  inline PhaseDataTpl<Scalar, Options, PhaseCollectionTpl>
  createData(const PhaseModelTpl<Scalar, Options, PhaseCollectionTpl> &phase_model)
  {
    return CreatePhaseData<Scalar, Options, PhaseCollectionTpl>::run(phase_model);
  }

} // namespace galileo

#endif // __galileo_core_phase_basic_visitors_hxx__
