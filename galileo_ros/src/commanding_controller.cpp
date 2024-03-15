#include "galileo_ros/commanding_controller.h"

// Constructor
CommandingController::CommandingController(ros::NodeHandle &node_handle, bool verbose) : node_handle_(node_handle), verbose_(verbose), udp_(UNITREE_LEGGED_SDK::LOWLEVEL, 8090, "192.168.12.1", 8007)
{
    // Initialize subscribers
    InitSubscribers();

    // Initialize publishers
    InitPublishers();

    // Initialize services
    InitServices();

    // Initialize WBC
    InitWBC();

    InitLowLevel();
}

void CommandingController::InitSubscribers()
{
}

void CommandingController::InitPublishers()
{
    // Publish to the Galileo command
    galileo_command_publisher_ = node_handle_.advertise<galileo_ros::RobotCommand>("legged_robot_command", 100);
}

void CommandingController::HardcodedPublisherInit()
{
    legged_robot_model_publisher_ = node_handle_.advertise<galileo_ros::ModelLocation>("legged_robot_model", 1);

    legged_parameter_location_publisher_ = node_handle_.advertise<std_msgs::String>("legged_parameter_location", 1);
}

void CommandingController::InitServices()
{
    // Get the desired state from the galileo solution
    get_desired_state_client_ = node_handle_.serviceClient<galileo_ros::DesiredStateInputCmd>("get_desired_state_input");

    // Subscribe to see if the problem is ready
    is_problem_ready_client_ = node_handle_.serviceClient<std_srvs::Trigger>("galileo_problem_is_ready");
}

void CommandingController::InitGalileo(galileo_ros::ModelLocation model_location_msg, std_msgs::String parameter_location_msg, galileo_ros::RobotCommand galileo_command)
{
    ros::Duration(0.3).sleep();

    // Call the is_problem_ready service
    while (!is_problem_ready_client_.exists() && ros::ok())
    {
        ros::Duration(0.1).sleep();
        ROS_INFO("Waiting for the galileo_problem_is_ready service to be ready...");
    }

    std_srvs::Trigger is_problem_ready_srv;
    is_problem_ready_client_.call(is_problem_ready_srv);
    while (!is_problem_ready_srv.response.success && ros::ok())
    {
        // Publish the messages
        legged_robot_model_publisher_.publish(model_location_msg);
        legged_parameter_location_publisher_.publish(parameter_location_msg);
        ros::Duration(0.1).sleep();
        is_problem_ready_client_.call(is_problem_ready_srv);
    }

    Eigen::VectorXd desired_state;
    Eigen::VectorXd desired_input;

    if (!getDesiredState(0, desired_state, desired_input))
    {
        // Publish the robot command
        galileo_command_publisher_.publish(galileo_command);
    }
}

void CommandingController::InitWBC()
{
    // For now we follow opem loop trajectory
}

void CommandingController::InitRobotState(const Eigen::VectorXd &desired_state)
{
    // publish robot state to gazebo
    AchieveLowLevelControl(desired_state, Eigen::VectorXd::Zero(desired_state.size()));
}

bool CommandingController::getDesiredState(float t_offset, Eigen::VectorXd &desired_state, Eigen::VectorXd &desired_input)
{
    // Create a new instance of the desired state input command
    galileo_ros::DesiredStateInputCmd desired_state_input_cmd;

    desired_state_input_cmd.request.time_offset_on_horizon = t_offset;

    // Call the service to get the desired state
    get_desired_state_client_.call(desired_state_input_cmd);

    // Check if the response is valid
    if (desired_state_input_cmd.response.valid_response)
    {
        // Get the desired state and input

        desired_state = Eigen::VectorXd::Map(desired_state_input_cmd.response.state_at_time_offset.data(), desired_state_input_cmd.response.state_at_time_offset.size());
        desired_input = Eigen::VectorXd::Map(desired_state_input_cmd.response.input_at_time_offset.data(), desired_state_input_cmd.response.input_at_time_offset.size());

        std::cout << "desired state gotten";
        return true;
    }
    else
    {
        return false;
    }
}

Eigen::VectorXd CommandingController::RunWBC(const Eigen::VectorXd &desired_state, const Eigen::VectorXd &desired_input, const Eigen::VectorXd &current_state) const
{
    // For now we follow opem loop trajectory
    return desired_input;
}

void CommandingController::RunNode()
{
    ros::Rate loop_rate(1000);
    Eigen::VectorXd desired_state;
    Eigen::VectorXd desired_input;

    getDesiredState(0, desired_state, desired_input);

    InitRobotState(desired_state);

    ros::Duration(0.3).sleep();

    std::time_t start_time = std::time(nullptr);
    while (ros::ok())
    {
        std::time_t current_time = std::time(nullptr);
        float t_offset = difftime(current_time, start_time);

        if (getDesiredState(t_offset, desired_state, desired_input))
        {
            AchieveLowLevelControl(desired_state, desired_input);
        }

        ros::spinOnce();
        loop_rate.sleep();
    }
}

void CommandingController::AchieveLowLevelControl(const Eigen::VectorXd &desired_state, const Eigen::VectorXd &desired_input)
{

    udp_.Recv();
    UNITREE_LEGGED_SDK::LowState state = {0};

    udp_.GetRecv(state);
    printf("%d  %f  %f\n", 0, state.motorState[0].q, state.motorState[0].dq);

    for (int joint_index = 0; joint_index < low_level_control_params_.leg_indeces.size(); joint_index++)
    {
        int motor_id = low_level_control_params_.leg_indeces[joint_index];
        int state_index = low_level_control_params_.state_indices[joint_index];
        cmd_.motorCmd[motor_id].q = 0; // desired_state(state_index);
        cmd_.motorCmd[motor_id].Kd = low_level_control_params_.kd;
        cmd_.motorCmd[motor_id].Kp = low_level_control_params_.kp;
        cmd_.motorCmd[motor_id].tau = low_level_control_params_.tau_offset[joint_index];
    }

    udp_.SetSend(cmd_);

    std::cout << "sending udp" << std::endl;
    udp_.Send();
}

int main(int argc, char **argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "go1_model_publisher_node");

    std::string model_file_location;
    if (argc > 1)
    {
        model_file_location = argv[1];
    }
    // Create a ROS node handle
    ros::NodeHandle node_handle;

    ros::Rate loop_rate(5);

    // Create a legged modelLocation message
    galileo_ros::ModelLocation model_location_msg;
    model_location_msg.model_file_location = model_file_location;
    model_location_msg.end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};
    model_location_msg.num_end_effectors = 4;

    // Create a parameter location message
    std_msgs::String parameter_location_msg;
    parameter_location_msg.data = "  "; // Empty right now

    std::vector<double> q0 = {
        0., 0., 0.339, 0., 0., 0., 1.,
        0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3};

    galileo_ros::RobotCommand robot_command_msg;
    robot_command_msg.initial_state = q0;
    robot_command_msg.target_state = q0;

    // Initialize ROS node
    ros::init(argc, argv, "commanding_controller_node");
    ros::NodeHandle nh;

    // Create an instance of CommandingController
    CommandingController commanding_controller(nh, true);

    commanding_controller.HardcodedPublisherInit();
    commanding_controller.InitGalileo(model_location_msg, parameter_location_msg, robot_command_msg);

    commanding_controller.RunNode();

    return 0;
}