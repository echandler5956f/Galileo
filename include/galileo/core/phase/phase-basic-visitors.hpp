#ifndef __galileo_core_phase_basic_visitors_hpp__
#define __galileo_core_phase_basic_visitors_hpp__

#include "galileo/core/phase/fwd.hpp"

namespace galileo
{

    template <typename Scalar, int Options, template <typename S, int O> class PhaseCollectionTpl>
    inline PhaseDataTpl<Scalar, Options, PhaseCollectionTpl>
    createData(const PhaseModelTpl<Scalar, Options, PhaseCollectionTpl> &phase_model);

} // namespace galileo

/* --- Details -------------------------------------------------------------------- */
// Included later
// #include "galileo/core/phase/core-basic-visitors.hxx"

#endif // __galileo_core_phase_basic_visitors_hpp__
