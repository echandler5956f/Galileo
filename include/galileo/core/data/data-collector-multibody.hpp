#ifndef __galileo_core_data_data_collector_multibody_hpp__
#define __galileo_core_data_data_collector_multibody_hpp__

#include <pinocchio/multibody/data.hpp>
#include "galileo/core/data/fwd.hpp"

namespace galileo
{

    // Pinocchio multibody data mixin
    template <typename Scalar, typename Derived>
    struct MultibodyDataMixin
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        pinocchio::DataTpl<Scalar> *pinocchio;

        MultibodyDataMixin(pinocchio::DataTpl<Scalar> *data) : pinocchio(data) {}
    };

} // namespace galileo

#endif // __galileo_core_data_data_collector_multibody_hpp__