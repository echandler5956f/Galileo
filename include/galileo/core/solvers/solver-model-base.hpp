#ifndef __galileo_core_solvers_solver_model_base_hpp__
#define __galileo_core_solvers_solver_model_base_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/core/phase/phase-generic.hpp"

namespace galileo
{

    /**
     * @brief Template class for representing a trajectory
     *
     * A trajectory is defined by a vector of segments, with continuity constraints
     * stitching each segment together. The list of segments can contain both
     * shooting and collocation types simultaneously.
     */
    template <typename _NumScalar, typename _VarScalar, template <typename, typename> class PhaseCollectionTpl>
    class TrajectoryModelTpl : CRTP<TrajectoryModelTpl<_NumScalar, _VarScalar, PhaseCollectionTpl>>
    {
    public:
    };

} // namespace galileo

#endif // __galileo_core_solvers_solver_model_base_hpp__