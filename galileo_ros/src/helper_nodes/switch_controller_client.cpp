// #include <ros/ros.h>
// #include "galileo_ros/DesiredStateInputCmd.h"

// int main(int argc, char **argv)
// {
//     // Initialize the ROS node
//     ros::init(argc, argv, "go1_poll_solution_node");

//     // Create a ROS node handle
//     ros::NodeHandle node_handle;

//     ros::ServiceClient client =
//         node_handle.serviceClient<galileo_ros::DesiredStateInputCmd>("get_desired_state_input");

//     galileo_ros::DesiredStateInputCmd srv;
//     // Request the horizon. For now it is 1.0

//     for (int i = 0; i < 50; i++)
//     {
//         srv.request.time_offset_on_horizon = i * 0.02;
//         std::cout << "calling service get_desired_state_input" << std::endl;
//         if (client.call(srv))
//         {
//             ROS_INFO("Desired state at time %f: ", srv.request.time_offset_on_horizon);
//             for (int i = 0; i < srv.response.state_at_time_offset.size(); i++)
//             {
//                 ROS_INFO("x[%d] = %f", i, srv.response.state_at_time_offset[i]);
//             }
//             ROS_INFO("Desired input at time %f: ", srv.request.time_offset_on_horizon);
//             for (int i = 0; i < srv.response.input_at_time_offset.size(); i++)
//             {
//                 ROS_INFO("u[%d] = %f", i, srv.response.input_at_time_offset[i]);
//             }
//         }
//         else
//         {
//             ROS_ERROR("Failed to call service get_desired_state_input");
//             std::cout << "failed" << std::endl;
//         }
//         ros::Duration(0.1).sleep();
//         std::cout << "called" << std::endl;
//         ros::spinOnce();
//     }
// }

#include <ros/ros.h>
#include <controller_manager_msgs/SwitchController.h>
#include <std_srvs/Empty.h>

int main(int argc, char **argv)
{
    ros::init(argc, argv, "switch_controller_client");
    ros::Duration(1.0).sleep();

    ros::NodeHandle nh;
    ros::ServiceClient client = nh.serviceClient<controller_manager_msgs::SwitchController>("/controller_manager/switch_controller");
    ros::ServiceClient unpause = nh.serviceClient<std_srvs::Empty>("/gazebo/unpause_physics");

    std_srvs::Empty empty;
    if (!unpause.call(empty)) {
        ROS_ERROR("Failed to unpause Gazebo");
        return 1;
    }

    controller_manager_msgs::SwitchController srv;
    srv.request.start_controllers.push_back("controllers/legged_controller");
    srv.request.stop_controllers.push_back("");
    srv.request.strictness = 0;
    srv.request.start_asap = false;
    srv.request.timeout = 0.0;

    if (client.call(srv))
    {
        ROS_INFO("Service call succeeded");
    }
    else
    {
        ROS_ERROR("Failed to call service");
        return 1;
    }

    return 0;
}