#include <ros/ros.h>
#include <std_msgs/String.h>

#include "galileo_ros/ModelLocation.h"
#include "galileo_ros/RobotSolution.h"
#include "galileo_ros/RobotCommand.h"

int main(int argc, char** argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "go1_model_publisher_node");

    std::string model_file_location;
    if(argc > 1){
        model_file_location = argv[1];
    }
    // Create a ROS node handle
    ros::NodeHandle node_handle;

    //Publish to legged_robot_model
    ros::Publisher legged_robot_model_publisher = 
        node_handle.advertise<galileo_ros::ModelLocation>("legged_robot_model", 1);

    //Publish to legged_parameter_location
    ros::Publisher legged_parameter_location_publisher = 
        node_handle.advertise<std_msgs::String>("legged_parameter_location", 1);

    ros::Publisher legged_robot_command_publisher = 
        node_handle.advertise<galileo_ros::RobotCommand>("legged_robot_command", 1);

    ros::Rate loop_rate(10);

    // Create a legged modelLocation message
    galileo_ros::ModelLocation model_location_msg;
    model_location_msg.model_file_location = model_file_location;
    model_location_msg.end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};
    model_location_msg.num_end_effectors = 4;

    // Create a parameter location message
    std_msgs::String parameter_location_msg;
    parameter_location_msg.data = "  "; // Empty right now

    std::vector<float> q0 = {
        0., 0., 0.339, 0., 0., 0., 1.,
        0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3};

    galileo_ros::RobotCommand robot_command_msg;
    robot_command_msg.initial_state = q0;
    robot_command_msg.target_state = q0;
    

    while(true){
        // Publish the messages
        legged_robot_model_publisher.publish(model_location_msg);
        legged_parameter_location_publisher.publish(parameter_location_msg);

        legged_robot_command_publisher.publish(robot_command_msg);

        // Spin and sleep
        ros::spinOnce();
    }



    return 0;
}