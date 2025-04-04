#ifndef __galileo_multibody_data_fwd_hpp__
#define __galileo_multibody_data_fwd_hpp__

#include <pinocchio/multibody/data.hpp>
#include "galileo/multibody/fwd.hpp"

namespace galileo
{

    namespace multibody
    {

        // Pinocchio multibody data mixin
        template <typename Scalar, typename Derived>
        struct MultibodyDataMixin
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            pinocchio::DataTpl<Scalar> *pinocchio;

            MultibodyDataMixin(pinocchio::DataTpl<Scalar> *data) : pinocchio(data) {}
        };

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_data_fwd_hpp__