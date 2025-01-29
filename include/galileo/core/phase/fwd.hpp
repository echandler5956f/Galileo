#ifndef __galileo_core_phase_fwd_hpp__
#define __galileo_core_phase_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    struct PhaseModelVoid
    {
    };

    struct PhaseDataVoid
    {
    };

    template <typename Scalar, int Options = context::Options>
    struct PhaseCollectionDefaultTpl;
    using PhaseCollectionDefault = PhaseCollectionDefaultTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
    struct PhaseModelTpl;
    using PhaseModel = PhaseModelTpl<context::Scalar>;

    template <
        typename Scalar,
        int Options = context::Options,
        template <typename S, int O> class PhaseCollectionTpl = PhaseCollectionDefaultTpl>
    struct PhaseDataTpl;
    using PhaseData = PhaseDataTpl<context::Scalar>;

} // namespace galileo

#include "galileo/core/fwd.hpp"

#endif // __galileo_core_phase_fwd_hpp__
