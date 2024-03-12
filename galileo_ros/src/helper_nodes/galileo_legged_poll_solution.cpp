#include <ros/ros.h>
#include "galileo_ros/DesiredStateInputCmd.h"

int main(int argc, char **argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "go1_poll_solution_node");

    // Create a ROS node handle
    ros::NodeHandle node_handle;

    ros::ServiceClient client =
        node_handle.serviceClient<galileo_ros::DesiredStateInputCmd>("get_desired_state_input");

    galileo_ros::DesiredStateInputCmd srv;
    // Request the horizon. For now it is 1.0

    for (int i = 0; i < 50; i++)
    {
        srv.request.time_offset_on_horizon = i * 0.02;
        std::cout << "calling service get_desired_state_input" << std::endl;
        if (client.call(srv))
        {
            ROS_INFO("Desired state at time %f: ", srv.request.time_offset_on_horizon);
            for (int i = 0; i < srv.response.state_at_time_offset.size(); i++)
            {
                ROS_INFO("x[%d] = %f", i, srv.response.state_at_time_offset[i]);
            }
            ROS_INFO("Desired input at time %f: ", srv.request.time_offset_on_horizon);
            for (int i = 0; i < srv.response.input_at_time_offset.size(); i++)
            {
                ROS_INFO("u[%d] = %f", i, srv.response.input_at_time_offset[i]);
            }
        }
        else
        {
            ROS_ERROR("Failed to call service get_desired_state_input");
            std::cout << "failed" << std::endl;
        }
        ros::Duration(0.1).sleep();
        std::cout << "called" << std::endl;
        ros::spinOnce();
    }
}
