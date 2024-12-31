#pragma once

#include "galileo/core/segment-base.hpp"

namespace galileo
{
    template <typename _Scalar>
    struct traits<TrajectoryModelTpl<_Scalar>>
    {
        typedef _Scalar Scalar;
    };

    /**
     * @brief Template class for representing a trajectory
     *
     * A trajectory is defined by a list of segments, with continuity constraints
     * stitching each segment together. The list of segments can contain both
     * shooting and collocation types simultaneously.
     */
    template <typename _Scalar>
    class TrajectoryModelTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
} // namespace galileo