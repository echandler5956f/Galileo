#ifndef __galileo_ros_legged_controller_hpp__
#define __galileo_ros_legged_controller_hpp__

#include <controller_interface/multi_interface_controller.h>
#include <hardware_interface/imu_sensor_interface.h>

#include <galileo/fwd.hpp>

#include "galileo_ros/common/hybrid-joint-interface.hpp"
#include "galileo_ros/common/contact-sensor-interface.hpp"

namespace galileo_ros
{

    class LeggedController : public controller_interface::MultiInterfaceController<HybridJointInterface,
                                                                                   hardware_interface::ImuSensorInterface,
                                                                                   ContactSensorInterface>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        bool init();
        void update(const ros::Time &time);
        void starting(const ros::Time &time);
        void stopping(const ros::Time &time);

        void setupStateEstimation();
        void updateStateEstimation(const ros::Time &time);

        void setupMPC();
        void updateMPC(const ros::Time &time);

        void setupWBC();
        void updateWBC(const ros::Time &time);

        void setupSimulator();
        void updateSimulator(const ros::Time &time);

        void setupVisualizer();
        void updateVisualizer(const ros::Time &time);

    }; // class LeggedController

} // namespace galileo_ros

#endif // __galileo_ros_legged_controller_hpp__