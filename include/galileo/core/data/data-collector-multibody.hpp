#ifndef __galileo_core_data_data_collector_multibody_hpp__
#define __galileo_core_data_data_collector_multibody_hpp__

namespace galileo
{

    // Pinocchio multibody data mixin
    template <typename Derived, typename PhaseSpec>
    struct MultibodyDataMixinTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using RobotData_t = typename PS::RobotData_t;

        RobotData_t *robot;

        MultibodyDataMixinTpl(RobotData_t *data)
            : robot(data) {}
    };

} // namespace galileo

#endif // __galileo_core_data_data_collector_multibody_hpp__
