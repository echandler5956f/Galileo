#ifndef __galileo_core_states_fwd_hpp__
#define __galileo_core_states_fwd_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar, typename NumScalar, int Options>
        class StateEuclidean;

        template <typename VarScalar, typename NumScalar, int Options>
        class StateSingleRigidBody;

        template <typename VarScalar, typename NumScalar, int Options>
        class StateCentroidalMomentumFullKinematics;
        
        template <typename VarScalar, typename NumScalar, int Options>
        class StateFloatingBase;

    } // namespace core

} // namespace galileo

#endif // __galileo_core_states_fwd_hpp__