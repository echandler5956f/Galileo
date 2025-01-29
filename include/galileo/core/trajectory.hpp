#pragma once

#include "galileo/core/segment/fwd.hpp"
#include "galileo/core/segment/segment-generic.hpp"

namespace galileo
{

    /**
     * @brief Template class for representing a trajectory
     *
     * A trajectory is defined by a vector of segments, with continuity constraints
     * stitching each segment together. The list of segments can contain both
     * shooting and collocation types simultaneously.
     */
    template <typename _Scalar, int _Options, template <typename, int> class SegmentCollectionTpl>
    class TrajectoryModelTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
} // namespace galileo