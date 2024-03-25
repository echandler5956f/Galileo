#include <ros/ros.h>
#include "galileo_ros/SolutionRequest.h"
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include <std_msgs/Float64.h>

std::vector<double> getSolutionAtTimeIdx(const galileo_ros::SolutionRequest &msg, int t_idx = 0)
{
    return std::vector<double>(msg.response.X_t_wrapped.begin() + msg.response.nx * t_idx, msg.response.X_t_wrapped.begin() + msg.response.nx * (1 + t_idx));
}

std::map<std::string, double> ConvertSolutionToJointMap(const galileo_ros::SolutionRequest &msg, int t_idx = 0)
{

    std::vector<double> Xt = getSolutionAtTimeIdx(msg, t_idx);

    std::vector<std::string> joint_names = msg.response.joint_names;

    std::vector<double> qj_t(Xt.begin() + msg.response.qj_index, Xt.begin() + msg.response.qj_index + joint_names.size());

    // map the joint names to the joint values
    std::map<std::string, double> joint_values;

    std::transform(joint_names.begin(), joint_names.end(), qj_t.begin(),
                   std::inserter(joint_values, joint_values.end()), std::make_pair<std::string const &, double const &>);

    return joint_values;
}

sensor_msgs::JointState getJointStateMessage(std::map<std::string, double> joint_values)
{
    sensor_msgs::JointState joint_state_msg;
    joint_state_msg.header.stamp = ros::Time::now();

    for (auto const &joint : joint_values)
    {
        joint_state_msg.name.push_back(joint.first);
        joint_state_msg.position.push_back(joint.second);
    }

    return joint_state_msg;
}

void setTransformationOn(const galileo_ros::SolutionRequest &msg, geometry_msgs::TransformStamped &body_transform)
{
    std::vector<double> transformation(msg.response.X_t_wrapped.begin() + msg.response.qj_index - 7, msg.response.X_t_wrapped.begin() + msg.response.qj_index);

    body_transform.header.stamp = ros::Time::now();
    body_transform.transform.translation.x = transformation[0];
    body_transform.transform.translation.y = transformation[1];
    body_transform.transform.translation.z = transformation[2];
    body_transform.transform.rotation.x = transformation[3];
    body_transform.transform.rotation.y = transformation[4];
    body_transform.transform.rotation.z = transformation[5];
    body_transform.transform.rotation.w = transformation[6];
}

int main(int argc, char **argv)
{
    double horizon = 0.6;
    // Initialize the ROS node
    ros::init(argc, argv, "galileo_ros_legged_rviz_node");
    ros::NodeHandle nh;
    std::string solver_id;
    if (!nh.getParam("galileo_ros/solver_id", solver_id))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    // Create a publisher to publish the solution
    ros::Publisher state_publisher = nh.advertise<sensor_msgs::JointState>("/joint_states", 1);

    ros::Publisher vis_time = nh.advertise<std_msgs::Float64>("/visualization_time", 1);

    // Create a ServiceClient to request the solution
    ros::ServiceClient solution_client = nh.serviceClient<galileo_ros::SolutionRequest>(solver_id + "_get_solution");
    geometry_msgs::TransformStamped body_transform;
    body_transform.header.frame_id = "world";
    body_transform.child_frame_id = "pelvis";

    // Create a TransformBroadcaster
    tf::TransformBroadcaster tf_broadcaster;
    std::cout << "Starting loop" << std::endl;

    // Main loop
    ros::Rate loop_rate(60);                 // 60 Hz
    ros::Time start_time = ros::Time::now(); // Get the start time
    while (ros::ok())
    {
        // Request the solution
        galileo_ros::SolutionRequest solution_request;

        ros::Time current_time = ros::Time::now();                 // Get the current time
        double elapsed_time = (current_time - start_time).toSec(); // Calculate the elapsed time in seconds
        elapsed_time = fmod(elapsed_time, horizon);                // TODO make this a parameter in the launch file
        solution_request.request.times = {elapsed_time};

        std_msgs::Float64 vis_time_msg;
        vis_time_msg.data = elapsed_time;
        vis_time.publish(vis_time_msg);

        if (solution_client.call(solution_request))
        {
            // std::cout << "Solution received" << std::endl;
            // Convert the solution to a map of joint values
            std::map<std::string, double> joint_values = ConvertSolutionToJointMap(solution_request);

            // get the horizon that the optimization is valid over
            horizon = solution_request.response.solution_horizon;

            // std::cout << "Publishing solution" << std::endl;

            // Create a JointState message
            sensor_msgs::JointState joint_state_msg = getJointStateMessage(joint_values);

            // std::cout << "Setting transformation" << std::endl;

            setTransformationOn(solution_request, body_transform);

            // std::cout << "Publishing transformation" << std::endl;

            // Publish the JointState message
            state_publisher.publish(joint_state_msg);

            // std::cout << "Getting transformation" << std::endl;

            // Broadcast the transform
            tf::StampedTransform tf_transform;
            tf::transformStampedMsgToTF(body_transform, tf_transform);

            // std::cout << "Broadcasting transformation" << std::endl;
            tf_broadcaster.sendTransform(tf_transform);
        }

        // Spin once to process callbacks and sleep to maintain the loop rate
        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}