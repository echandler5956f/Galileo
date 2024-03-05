#include <ros/ros.h>

#include "galileo_ros/ModelLocation.h"
#include "galileo_ros/RobotSolution.h"

int main(int argc, char** argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "go1_model_publisher_node");

    // Create a ROS node handle
    ros::NodeHandle node_handle;

    //Publish to legged_robot_model
    ros::Publisher legged_robot_model_publisher = 
        node_handle.advertise<galileo_ros::ModelLocation>("legged_robot_model", 1);

    //Publish to legged_parameter_location
    ros::Publisher legged_parameter_location_publisher = 
        node_handle.advertise<std_msgs::String>("legged_parameter_location", 1);

    // Create a legged modelLocation message
    galileo_ros::ModelLocation model_location_msg;
    model_location_msg.model_location = "/home/akshay/Documents/BiQu/src/Galileo/resources/go1/urdf/go1.urdf";
    model_location_msg.end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};

    // Create a parameter location message
    std_msgs::String parameter_location_msg;
    parameter_location_msg.data = "  "; // Empty right now

    while(true){
        // Publish the messages
        legged_robot_model_publisher.publish(model_location_msg);
        legged_parameter_location_publisher.publish(parameter_location_msg);
        
        ros::spinOnce();
    }



    return 0;
}