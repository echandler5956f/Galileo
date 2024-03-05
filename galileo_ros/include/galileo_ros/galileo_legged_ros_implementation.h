

#pragma once

#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/Twist.h>
#include <ros/subscriber.h>

#include <Eigen/Core>

#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"
#include "galileo/opt/TrajectoryOpt.h"

#include "galileo_ros/ModelLocation.h"
#include "galileo_ros/RobotSolution.h"


#define BASE_FRAME "base"


/**
 * @class GalileoLeggedROSImplementation
 * @brief This class represents the implementation of the ROS interface for the Galileo Legged robot.
 * 
 * The GalileoLeggedROSImplementation class handles calls to Galileo.
 */
class GalileoLeggedROSImplementation {
    public: 

    // casadi DM of size nx
    typedef T_ROBOT_STATE casadi::DM;
    
    using LeggedConstraintBuilderType = std::shared_ptr<gallileo::legged::constraints::ConstraintBuilder<gallileo::legged::constraints::LeggedRobotProblemData>>; 
    using LeggedTrajOpt = galileo::opt::TrajectoryOpt<galileo::legged::LeggedRobotProblemData, galileo::legged::contact::ContactMode>;

    
    /**
     * @brief Construct a new Galileo Legged ROS Implementation object to handle calls to Galileo.
     * 
     * @param node_handle The ROS node handle.
     */
    GalileoLeggedROSImplementation(::ros::NodeHandle& node_handle);

    /**
     * @brief Loads the robot model from the given model file and sets the end effectors.
     * 
     * @param model_file The path to the robot model file.
     * @param end_effector_names The names of the end effectors.
     */
    void LoadModelCallback(const galileo_ros::ModelLocation::ConstPtr& msg);

    /**
     * @brief Updates the parameters from the given parameter file, and create the costs.
     * 
     * @param parameter_file The path to the parameter file.
     */
    void ParameterCallback(const std_msgs::String::ConstPtr& msg);

    /**
     * @brief Publishes the last solution.
     */
    void PublishLastSolution();
-
    /**
     * @brief Updates the solution based on the initial state.
     * 
     * @param initial_state The initial state of the robot.
     */
    void UpdateSolution( T_ROBOT_STATE initial_state );    
    
    
    /**
     * @brief Creates the problem data from the loaded model/parameters.
    */
    void CreateProblemData();

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
    void InitServices(){} // no services yet.
        
    /**
     * @brief Loads the robot model from the given model file and sets the end effectors.
     * 
     * @param model_file The path to the robot model file.
     * @param end_effector_names The names of the end effectors.
     */
    void LoadModel(const std::string& model_file, const std::vector<std::string>& end_effector_names);

    // Hardcoded functions for the Go1 implementation. To be changed.

    /**
     * @brief Loads the parameters from the given parameter file.
     * 
     * @param parameter_file The path to the parameter file.
     */
    void LoadParameters(const std::string& parameter_file);
    


    /**
     * @brief Updates the robot running and terminal costs L and Phi. 
     * 
    */
    void CreateCost( casadi::Function &L, casadi::Function &Phi ) const;


    
    /**
     * @brief Gets the legged constraint builders.
    */
    std::vector< LeggedConstraintBuilderType > 
        getLeggedConstraintBuilders() const;



    std::shared_ptr<galileo::legged::LeggedBody> robot_; /**< The robot model. */

    std::shared_ptr<galileo::legged::LeggedRobotStates> states_; /**< Definition of the state. */

    std::shared_ptr<galileo::legged::LeggedRobotProblemData> problem_data_; /**< The problem data. */

    std::shared_ptr<LeggedTrajOpt> trajectory_opt_; /**< The trajectory optimizer. */

    std::shared_ptr<galileo::environment::EnvironmentSurfaces> surfaces_; /**< The surfaces. */
    
    std::shared_ptr<galileo::legged::contact::ContactSequence> contact_sequence_; /**< The gait. */



    ::ros::NodeHandle& node_handle_; /**< The ROS node handle. */

    ::ros::Subscriber initial_state_subscriber_; /**< The subscriber for the initial state. */

    ::ros::Subscriber target_pose_subscriber_; /**< The subscriber for the target pose. */

    ::ros::Publisher solution_publisher_; /**< The publisher for the solution. */

    ::ros::Subscriber parameter_location_subscriber_; /**< The subscriber for the parameter file location. */

    ::ros::Subscriber robot_model_subscriber_; /**< The subscriber for the robot model. */



};