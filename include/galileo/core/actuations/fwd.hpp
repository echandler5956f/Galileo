#ifndef __galileo_core_actuations_fwd_hpp__
#define __galileo_core_actuations_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationModelFullTpl;
        template <typename VarScalar, typename NumScalar, int Options>
        struct ActuationDataFullTpl;

        template <typename VarScalar, typename NumScalar, int Options>
        class ActuationModelFloatingBaseTpl;
        template <typename VarScalar, typename NumScalar, int Options>
        struct ActuationDataFloatingBaseTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_fwd_hpp__