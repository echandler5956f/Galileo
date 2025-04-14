#ifndef __galileo_core_data_data_collector_multibody_hpp__
#define __galileo_core_data_data_collector_multibody_hpp__

#include "galileo/core/data/fwd.hpp"

namespace galileo
{

    // Pinocchio multibody data mixin
    template <typename Derived, typename PhaseSpec>
    struct MultibodyDataMixinTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        typename PS::RobotData_t *robot;

        MultibodyDataMixinTpl(typename PS::RobotData_t *data) : robot(data) {}
        
    };

} // namespace galileo

#endif // __galileo_core_data_data_collector_multibody_hpp__