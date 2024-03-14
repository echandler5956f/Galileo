#include <ros/ros.h>
#include "galileo_ros/DesiredStateInputCmd.h"
#include "sensor_msgs/JointState.h"

int main(int argc, char **argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "go1_poll_solution_node");

    sleep(2);

    // Create a ROS node handle
    ros::NodeHandle node_handle;

    ros::ServiceClient client =
        node_handle.serviceClient<galileo_ros::DesiredStateInputCmd>("get_desired_state_input");

    // Publisher to publish jointstate msg for Rviz visual
    ros::Publisher desiredJointState_pub =
        node_handle.advertise<sensor_msgs::JointState>("joint_states", 100);

    galileo_ros::DesiredStateInputCmd srv;

    // Request the horizon. For now it is 1.0

    // q = [hx, hy, hz, lx, ly, lz,
    //      qb_pos_x, qb_pos_y, qb_pos_z, qb_ang_x, qb_ang_y, qb_ang_z,
    //      qj0... qj11]
    // ["FL_hip_joint", "FL_thigh_joint", "FL_calf_joint", "FR_hip_joint", "FR_thigh_joint", "FR_calf_joint", "RL_hip_joint", "RL_thigh_joint", "RL_calf_joint", "RR_hip_joint", "RR_thigh_joint", "RR_calf_joint"]

    sensor_msgs::JointState desiredJS_msg;
    desiredJS_msg.name = {"FL_hip_joint", "FL_thigh_joint", "FL_calf_joint", "FR_hip_joint", "FR_thigh_joint", "FR_calf_joint", "RL_hip_joint", "RL_thigh_joint", "RL_calf_joint", "RR_hip_joint", "RR_thigh_joint", "RR_calf_joint"};
    std::cout << "calling service get_desired_state_input" << std::endl;
    for (int i = 0; i < 100; i++)
    {
        desiredJS_msg.position.clear();
        srv.request.time_offset_on_horizon = i * 0.01;

        if (client.call(srv))
        {
            // ROS_INFO("Desired state at time %f: ", srv.request.time_offset_on_horizon);
            for (int i = 0; i < srv.response.state_at_time_offset.size(); i++)
            {
                if (i >= 12)
                {
                    desiredJS_msg.position.push_back(srv.response.state_at_time_offset[i]);
                    ROS_INFO("x[%d] = %f", i, srv.response.state_at_time_offset[i]);
                }
            }
            // ROS_INFO("Desired input at time %f: ", srv.request.time_offset_on_horizon);
            // for (int i = 0; i < srv.response.input_at_time_offset.size(); i++)
            // {
            //     ROS_INFO("u[%d] = %f", i, srv.response.input_at_time_offset[i]);
            // }
        }
        else
        {
            ROS_ERROR("Failed to call service get_desired_state_input");
            std::cout << "failed" << std::endl;
        }
        // ROS_INFO("Sending %d Joint States to /joint_states:", desiredJS_msg.position.size());

        desiredJS_msg.header.stamp = ros::Time::now();
        desiredJointState_pub.publish(desiredJS_msg);

        ros::Duration(0.1).sleep();
        ros::spinOnce();
    }
    std::cout << "Finished publishing" << std::endl;
    return 0;
}
