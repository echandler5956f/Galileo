#include "galileo_ros/GalileoLeggedRos.h"

#include <ros/ros.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "galileo_ros_legged_solver_node");

    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();
    std::string solver_id;

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    ROS_INFO("Creating the solver object\n");
    galileo::legged::GalileoLeggedRos galileo_legged_solver(nh, solver_id);

    ros::spin();

    return 0;
}