#ifndef __galileo_core_actuations_fwd_hpp__
#define __galileo_core_actuations_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationModelFull;

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationDataFull;

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationModelFloatingBase;

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationDataFloatingBase;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_fwd_hpp__