#pragma once

#include <pinocchio/fwd.hpp>

#include <ros/subscriber.h>
#include <ros/publisher.h>
#include <ros/service.h>
// Include the necessary ROS header files
#include <ros/ros.h>
#include <std_msgs/String.h> // Assuming the parameter location strings are published as std_msgs::String
#include <std_srvs/Trigger.h>

#include <Eigen/Core>

#include <thread>

#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedRobotStates.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"
#include "galileo/opt/TrajectoryOpt.h"

#include "galileo_ros/ModelLocation.h"
#include "galileo_ros/RobotCommand.h"
#include "galileo_ros/DesiredStateInputCmd.h"

class CommandingController
{
public:
    CommandingController(::ros::NodeHandle &node_handle, bool verbose = false);

    /**
     * @brief Publish specifically for LeggedROS Go1
     */
    void HardcodedPublisherInit();

    /**
     * @brief Initializes the Galileo Legged ROS Implementation and attempts to create a solution.
     */
    void InitGalileo(galileo_ros::ModelLocation model_location_msg, std_msgs::String parameter_location_msg, galileo_ros::RobotCommand galileo_command);

    /**
     *
     */
    void RunNode();

    /**
     * @brief Create and publish the robot command to the robot. This should be threaded.
     *
     */
    void PublishRobotJointCommand(float t_offset);

private:
    /**
     * @brief Initializes the subscribers for the Galileo Legged ROS Implementation.
     */
    void InitSubscribers();

    /**
     * @brief Initializes the publishers for the Galileo Legged ROS Implementation.
     */
    void InitPublishers();

    /**
     * @brief Initializes the services for the Galileo Legged ROS Implementation.
     */
    void InitServices();

    /**
     * @brief Initializes the whole body controller for the Galileo Legged ROS Implementation.
     */
    void InitWBC();

    /**
     * @brief Initializes the robot state in gazebo before running the WBC to follow the open-loop trajectory.
     */
    void InitRobotState();

    /**
     * @brief get the desired state at time t from the galileo solution
     */
    bool getDesiredState(float t_offset, Eigen::VectorXd &desired_state, Eigen::VectorXd &desired_input);

    /**
     * @brief get the desired tau from the WBC
     */
    Eigen::VectorXd RunWBC(const Eigen::VectorXd &desired_state,
                           const Eigen::VectorXd &desired_input,
                           const Eigen::VectorXd &current_state) const;

    ::ros::NodeHandle &node_handle_; /**< The ROS node handle. */

    ::ros::ServiceClient get_desired_state_client_; /**< The ROS service client for getting the desired state. */

    ::ros::ServiceClient is_problem_ready_client_; /**< The ROS service client for checking if the problem is ready. */

    ::ros::Publisher galileo_command_publisher_; /**< The ROS publisher for the Galileo command. */

    ::ros::Publisher robot_state_publisher_; /**< The ROS publisher for the robot state. */

    ::ros::Publisher tau_publisher_; /**< The ROS publisher for the desired tau. */

    ::ros::Subscriber robot_state_subscriber_; /**< The ROS subscriber for the robot state. */

    // Publish to legged_robot_model
    ros::Publisher legged_robot_model_publisher_;

    // Publish to legged_parameter_location
    ros::Publisher legged_parameter_location_publisher_;

    // std::shared_ptr<WbcBase> wbc_; /**< The whole body controller. */

    bool verbose_; /**< Whether to print verbose output. */

    struct WBCParameters
    {
        std::vector<int> joint_indices;
        std::vector<float> max_tau;
        std::vector<float> min_tau;

        float kp;
    };
};