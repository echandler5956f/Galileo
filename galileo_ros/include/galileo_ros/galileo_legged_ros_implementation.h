

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

#include "galileo/math/Quat2Euler.h"
#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedRobotStates.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"
#include "galileo/opt/TrajectoryOpt.h"

#include "galileo_ros/ModelLocation.h"
#include "galileo_ros/RobotSolution.h"
#include "galileo_ros/RobotCommand.h"
#include "galileo_ros/DesiredStateInputCmd.h"

#define BASE_FRAME "base"

/**
 * @class GalileoLeggedROSImplementation
 * @brief This class represents the implementation of the ROS interface for the Galileo Legged robot.
 *
 * The GalileoLeggedROSImplementation class handles calls to Galileo.
 */
class GalileoLeggedROSImplementation
{
public:
    // casadi DM of size nx
    using T_ROBOT_STATE = casadi::DM;

    using LeggedRobotProblemData = galileo::legged::constraints::LeggedRobotProblemData;
    using GeneralProblemData = galileo::opt::GeneralProblemData;

    using LeggedConstraintBuilderType = std::shared_ptr<galileo::opt::ConstraintBuilder<LeggedRobotProblemData>>;
    using LeggedTrajOpt = galileo::opt::TrajectoryOpt<LeggedRobotProblemData, galileo::legged::contact::ContactMode>;

    /**
     * @brief Construct a new Galileo Legged ROS Implementation object to handle calls to Galileo.
     *
     * @param node_handle The ROS node handle.
     */
    GalileoLeggedROSImplementation(::ros::NodeHandle &node_handle, bool verbose = false);

    /**LeggedConstraintBuilderType
     * @brief Loads the robot model from the given model file and sets the end effectors.
     *
     * @param model_file The path to the robot model file.
     * @param end_effector_names The names of the end effectors.
     */
    void LoadModelCallback(const galileo_ros::ModelLocation::ConstPtr &msg);

    /**
     * @brief Updates the parameters from the given parameter file, and create the costs.
     *
     * @param parameter_file The path to the parameter file.
     */
    void ParameterCallback(const std_msgs::String::ConstPtr &msg);

    /**
     * @brief Updates the solution with a new initial state.
     * Currently DOES NOT update a target state.
     *
     * @param msg The initial state of the robot.
     */
    void RobotCommandCallback(const galileo_ros::RobotCommand::ConstPtr &msg);

    /**
     * @brief returns the desired state and input at time "t"
     *
     */
    bool DesiredStateInputCallback(galileo_ros::DesiredStateInputCmd::Request &req,
                                   galileo_ros::DesiredStateInputCmd::Response &res);

    bool IsReadyCallback(std_srvs::Trigger::Request &req,
                         std_srvs::Trigger::Response &res)
    {
        res.success = (bool)trajectory_opt_;
        return true;
    }

    /**
     * @brief Publishes the last solution.
     */
    void PublishLastSolution();
    /**
     * @brief Updates the solution based on the initial state.
     *
     * @param initial_state The initial state of the robot.
     */
    void UpdateSolution(T_ROBOT_STATE initial_state);

    /**
     * @brief Creates the problem data from the loaded model/parameters.
     */
    void CreateProblemData(T_ROBOT_STATE initial_state);

    /**
     * @brief Creates the trajectory optimizer if the problem data is set.
     */
    void CreateTrajOptSolver();

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
     * @brief Loads the robot model from the given model file and sets the end effectors.
     *
     * @param model_file The path to the robot model file.
     * @param end_effector_names The names of the end effectors.
     */
    void LoadModel(const std::string &model_file, const std::vector<std::string> &end_effector_names);

    // Hardcoded functions for the Go1 implementation. To be changed.

    /**
     * @brief Loads the parameters from the given parameter file.
     *
     * @param parameter_file The path to the parameter file.
     */
    void LoadParameters(const std::string &parameter_file);

    /**
     * @brief Updates the robot running and terminal costs L and Phi.
     *
     */
    void CreateCost(casadi::Function &L, casadi::Function &Phi) const;

    /**
     * @brief Part of the hardcoded initialization. Returns X0.
     */
    casadi::DM getX0(double *q0) const;

    /**
     * @brief Gets the legged constraint builders.
     */
    std::vector<LeggedConstraintBuilderType>
    getLeggedConstraintBuilders() const;

    std::shared_ptr<galileo::legged::LeggedBody> robot_; /**< The robot model. */

    std::shared_ptr<galileo::legged::LeggedRobotStates> states_; /**< Definition of the state. */

    std::shared_ptr<LeggedRobotProblemData> problem_data_; /**< The problem data. */

    std::shared_ptr<LeggedTrajOpt> trajectory_opt_; /**< The trajectory optimizer. */

    std::shared_ptr<galileo::legged::environment::EnvironmentSurfaces> surfaces_; /**< The surfaces. */

    std::shared_ptr<galileo::legged::contact::ContactSequence> contact_sequence_; /**< The gait. */

    casadi::Dict opts_;

    std::shared_ptr<galileo::opt::DecisionDataBuilder<LeggedRobotProblemData>> decision_builder_;

    bool verbose_; /**< Whether to print verbose output. */

    ::ros::NodeHandle &node_handle_; /**< The ROS node handle. */

    ::ros::Subscriber initial_state_subscriber_; /**< The subscriber for the initial state. */

    ::ros::Subscriber target_pose_subscriber_; /**< The subscriber for the target pose. */

    ::ros::Subscriber parameter_location_subscriber_; /**< The subscriber for the parameter file location. */

    ::ros::Subscriber robot_model_subscriber_; /**< The subscriber for the robot model. */

    std::thread galileo_solver_; /**< The thread for running the optimizer */

    std::mutex solution_mutex_; /**< Bar the solution from being written to/accessed. */

    std::mutex optimizing_mutex_; /**< Bar the solution from being written to/accessed. */

    /** For getting solutions. We eventually want to publish the last known solution at a high frequency */
    ::ros::Publisher solution_publisher_; /**< Publishes the last solution. */

    ::ros::Subscriber robot_command_subscriber_; /**< Updates the solution when called. */

    ::ros::ServiceServer desired_state_input_service_; /**< Returns the desired state and input at time "t". */

    ::ros::ServiceServer problem_ready_; /**< Returns whether the problem is ready to be "updated". */

    casadi::MXVector solution_; /**< The last solution. */

    bool fully_initted_ = false; /**< Whether a solution has been made. */
};