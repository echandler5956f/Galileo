#pragma once

#include "galileo/core/segment-base.hpp"

namespace galileo
{
    template <typename Derived, typename Scalar>
    class TrajectoryModelAbstractTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef MathBaseTpl<Scalar> MathBase;
        typedef StateAbstractTpl<Scalar> StateAbstract;
        typedef typename MathBase::VectorXs VectorXs;
        typedef typename MathBase::MatrixXs MatrixXs;
    };
}