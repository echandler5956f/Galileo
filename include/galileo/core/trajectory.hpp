#ifndef __galileo_core_trajectory_hpp__
#define __galileo_core_trajectory_hpp__

#include "galileo/core/fwd.hpp"
#include "galileo/core/phase/phase-generic.hpp"

namespace galileo
{

    template <typename _NumScalar, typename _VarScalar, template <typename, typename> class PhaseCollectionTpl>
    struct traits<TrajectoryModelTpl<_NumScalar, _VarScalar, PhaseCollectionTpl>>
    {
        using NumScalar = _NumScalar;
        using VarScalar = _VarScalar;
    };

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
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
        using NumScalar = _NumScalar;
        using VarScalar = _VarScalar;
    
        using PhaseModel = PhaseModelTpl<NumScalar, VarScalar, PhaseCollectionTpl>;
        using PhaseData = PhaseDataTpl<NumScalar, VarScalar, PhaseCollectionTpl>;

        using PhaseModelVector = GALILEO_ALIGNED_STD_VECTOR(PhaseModel);
        using PhaseDataVector = GALILEO_ALIGNED_STD_VECTOR(PhaseData);

        PhaseModelVector phases_;
    };

} // namespace galileo

#endif // __galileo_core_trajectory_hpp__