#include "galileo_ros/GalileoLeggedRos.h"

#include <ros/ros.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "galileo_ros_legged_solver_node");

    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();
    std::string solver_id;
    std::string sol_data_dir;
    std::string plot_dir;

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/sol_data_dir", sol_data_dir) ||
        !nh->getParam("galileo_ros/plot_dir", plot_dir))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    ROS_INFO("Creating the solver object\n");
    galileo::legged::GalileoLeggedRos galileo_legged_solver(nh, solver_id, sol_data_dir, plot_dir);

    ros::spin();

    return 0;
}