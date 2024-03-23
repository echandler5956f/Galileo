#include "galileo_ros/GalileoLeggedRos.h"

#include <ros/ros.h>

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        std::cerr << "Usage: rosrun galileo_ros galileo_legged_test_node <solver_id> " << std::endl;
        return 1;
    }

    std::string solver_id = argv[1];

    ros::init(argc, argv, solver_id + "_solver_node");
    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    galileo::legged::GalileoLeggedRos galileo_legged_solver(nh, solver_id);

    ros::spin();

    return 0;
}