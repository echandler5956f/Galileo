#ifndef __galileo_core_costs_fwd_hpp__
#define __galileo_core_costs_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, bool ShareData = false, int NX = -1, int NU = -1, int Options>
        struct CostMatricesTpl; // forward declaration

        template <typename VarScalar, typename NumScalar, int Options>
        using CostMatricesDynamic = CostMatricesTpl<VarScalar, NumScalar, Options>;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_costs_fwd_hpp__