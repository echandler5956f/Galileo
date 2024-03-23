#include <ros/ros.h>
#include "galileo_ros/SolutionRequest.h"
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>

std::map<std::string, double> ConvertSolutionToJointMap(const galileo_ros::SolutionRequest &msg)
{
    // For now we only concern ourself with the first state
    std::vector<double> Xt = msg.X_t_wrapped(msg.X_t_wrapped.begin(), msg.X_t_wrapped.begin() + msg.nx);

    std::vector<std::string> joint_names = msg.joint_names;

    std::vector<double> qj_t = Xt(Xt.begin() + msg.qj_index, Xt.begin() + msg.qj_index + joint_names.size());

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
    std::vector<double> transformation = msg.X_t_wrapped(msg.X_t_wrapped.begin() + msg.q_index, msg.X_t_wrapped.begin() + msg.q_index + 7);

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
    // Initialize the ROS node
    ros::init(argc, argv, "galileo_legged_rvis_node");
    ros::NodeHandle nh;

    if (argc != 2)
    {
        ROS_ERROR("Usage: rosrun galileo_ros galileo_legged_rvis_node <solver_id>");
        return 1;
    }
    std::string solver_id = argv[1];

    // Create a publisher to publish the solution
    //
    ros::Publisher state_publisher = nh.advertise<sensor_msgs::JointState>("/joint_states", 1);
    ros::Publisher tf_broadcaster = nh.advertise<tf::TransformBroadcaster>("/tf", 1);

    // Create a ServiceClient to request the solution

    ros::ServiceClient solution_client = nh.serviceClient<galileo_ros::SolutionRequest>(solver_id + "_get_solution");
    geometry_msgs::TransformStamped body_transform;
    body_transform.header.frame_id = "odom";
    body_transform.child_frame_id = "body";

    // Main loop
    ros::Rate loop_rate(50); // 50 Hz
    while (ros::ok())
    {
        // Request the solution
        galileo_ros::SolutionRequest solution_request;

        // SET TO 0 FOR NOW
        solution_request.request.times = {0};

        if (solution_client.call(solution_request))
        {
            // Convert the solution to a map of joint values
            std::map<std::string, double> joint_values = ConvertSolutionToJointMap(solution_request.response);

            // Create a JointState message
            sensor_msgs::JointState joint_state_msg = getJointStateMessage(joint_values);

            setTransformationOn(solution_request.response, body_transform);

            // Publish the JointState message
            state_publisher.publish(joint_state_msg);
            tf_broadcaster.sendTransform(body_transform);
        }

        // Spin once to process callbacks and sleep to maintain the loop rate
        ros::spinOnce();
        loop_rate.sleep();
    }

    return 0;
}