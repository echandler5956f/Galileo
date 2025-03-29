#ifndef __galileo_core_states_fwd_hpp__
#define __galileo_core_states_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, int Options, int NX = -1, int NU = -1, int NDX = -1>
        class StateEuclideanTpl;

        template <typename VarScalar, typename NumScalar, int Options, int NQ = -1, int NV = -1, int NFb = -1>
        class StateMultibodyTpl;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_fwd_hpp__