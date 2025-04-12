#ifndef __galileo_predictive_mrt_base_hpp__
#define __galileo_predictive_mrt_base_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/trajectory.hpp"

#define GALILEO_MRT_BASIC_TYPEDEF(MRT)                 \
    using VarScalar = typename traits<MRT>::VarScalar; \
    using NumScalar = typename traits<MRT>::NumScalar; \
    static constexpr int Options = traits<MRT>::Options;

namespace galileo
{

    template <class Derived>
    class MRTBase : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using MRTBaseDerived = typename traits<Derived>::MRTBaseDerived;
        GALILEO_MRT_BASIC_TYPEDEF(MRTBaseDerived);

    }; // class MRTBase

} // namespace galileo

#endif // __galileo_predictive_mrt_base_hpp__