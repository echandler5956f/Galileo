#include "galileo_ros/GalileoLeggedRos.h"

std::vector<double> getX0(int nx, int q_index)
{
    galileo::legged::ConfigVector q0_vec = (galileo::legged::ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

    std::vector<double> X0;
    for (int j = 0; j < nx; j++)
    {
        X0.push_back(0);
    }

    for (int j = 0; j < 19; j++)
    {
        X0[q_index + j] = q0_vec(j);
    }

    return X0;
}

std::vector<double> getXf(int nx, int q_index)
{
    std::vector<double> Xf = getX0(nx, q_index);
    Xf[q_index] = 0.2;
    return Xf;
}

void getProblemDataMessages(std::string urdf_name, std::string parameter_file_name,
                            galileo_ros::RobotModel &robot_model_cmd,
                            galileo_ros::ParameterFileLocation &parameter_location_cmd,
                            galileo_ros::ContactSequence &contact_sequence_cmd, galileo_ros::GalileoCommand &galileo_cmd_msg)
{
    std::string model_location = urdf_name;
    std::string parameter_location = parameter_file_name;

    robot_model_cmd.model_file_location = model_location;

    // std::string end_effectors[] = {"FL_foot", "FR_foot", "RL_foot", "RR_foot"};
    std::string end_effectors[] = {"l_foot_v_ft_link", "r_foot_v_ft_link"};
    robot_model_cmd.end_effector_names = std::vector<std::string>(end_effectors, end_effectors + sizeof(end_effectors) / sizeof(std::string));

    parameter_location_cmd.parameter_file_location = parameter_location;

    std::vector<int> knot_num = {25, 16, 8, 16, 25};
    std::vector<double> knot_time = {0.25, 0.1667, 0.08, 0.1667, 0.25};
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
        {0, 0},
        {-1, 0},
        {0, 0},
        {0, -1},
        {0, 0}}; // static walk

    // std::vector<int> knot_num = {5, 30, 30, 30, 30, 5};
    // std::vector<double> knot_time = {0.05, 0.3, 0.3, 0.3, 0.3, 0.05};
    // std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
    //     {0, 0, 0, 0},
    //     {0, -1, -1, 0},
    //     {-1, 0, 0, -1},
    //     {0, -1, -1, 0},
    //     {-1, 0, 0, -1},
    //     {0, 0, 0, 0}}; // trot

    // std::vector<int> knot_num = {30, 30, 30};
    // std::vector<double> knot_time = {0.15, 0.3, 0.15};
    // std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
    //     {0, 0, 0, 0},
    //     {-1, -1, -1, -1},
    //     {0, 0, 0, 0}}; // jump/

    // std::vector<int> knot_num = {50};
    // std::vector<double> knot_time = {0.6};
    // std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
    //     {0, 0, 0, 0}}; // stand

    for (int phase_number = 0; phase_number < knot_num.size(); phase_number++)
    {
        galileo_ros::ContactPhase phase;
        phase.knot_num = knot_num[phase_number];
        phase.knot_time = knot_time[phase_number];
        phase.contact_surface_ids = contact_surfaces[phase_number];
        contact_sequence_cmd.phases.push_back(phase);
    }

    galileo::legged::LeggedBody robot(model_location, robot_model_cmd.end_effector_names.size(), end_effectors);

    std::vector<double> X0 = getX0(robot.si->nx, robot.si->q_index);
    std::vector<double> Xf = getXf(robot.si->nx, robot.si->q_index);

    galileo_cmd_msg.command_type = "init";

    galileo_cmd_msg.X_initial = X0;
    galileo_cmd_msg.X_goal = Xf;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "galileo_ros_legged_initialization_node");
    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    std::string solver_id;
    std::string urdf_file_name;
    std::string parameter_file_name;

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/urdf_filename", urdf_file_name) ||
        !nh->getParam("galileo_ros/parameters_location", parameter_file_name))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    galileo_ros::RobotModel model_location_msg;
    galileo_ros::ParameterFileLocation parameter_location_msg;
    galileo_ros::ContactSequence contact_sequence_msg;
    galileo_ros::GalileoCommand galileo_cmd_msg;
    galileo_ros::EnvironmentSurface surface_msg;

    ROS_INFO("Generating problem data messages\n");

    getProblemDataMessages(urdf_file_name, parameter_file_name,
                           model_location_msg,
                           parameter_location_msg,
                           contact_sequence_msg, galileo_cmd_msg);

    ROS_INFO("Publishing model location, parameter location and contact sequence\n");

    ros::Publisher model_location_pub = nh->advertise<galileo_ros::RobotModel>(solver_id + "_model_location", 1);
    ros::Publisher parameter_location_pub = nh->advertise<galileo_ros::ParameterFileLocation>(solver_id + "_parameter_location", 1);
    ros::Publisher contact_sequence_pub = nh->advertise<galileo_ros::ContactSequence>(solver_id + "_contact_sequence", 1);
    ros::Publisher command_publisher = nh->advertise<galileo_ros::GalileoCommand>(solver_id + "_command", 1);
    ros::Publisher surface_pub = nh->advertise<galileo_ros::EnvironmentSurface>(solver_id + "_add_environment_surface", 1);
    ros::ServiceClient init_state_client = nh->serviceClient<galileo_ros::InitState>(solver_id + "_init_state_service");

    galileo_ros::InitState init_state;
    init_state.request.call = true;

    ROS_INFO("Calling init state\n");
    ros::spinOnce();
    while (!init_state_client.call(init_state))
    {
        ROS_INFO("Still initing state\n");
        ros::spinOnce();
        ros::Duration(0.3).sleep();
    }

    ROS_INFO("Called init state\n");

    while (true)
    {
        init_state_client.call(init_state);
        ROS_INFO("Model set: %d, solver parameters set: %d, contact sequence set: %d, fully initted: %d\n",
                 init_state.response.model_set,
                 init_state.response.solver_parameters_set,
                 init_state.response.contact_sequence_set,
                 init_state.response.fully_initted);

        if (!init_state.response.model_set)
        {

            ROS_INFO("Publishing model location\n");
            model_location_pub.publish(model_location_msg);
        }
        else if (!init_state.response.solver_parameters_set)
        {
            ROS_INFO("Publishing parameter location\n");
            parameter_location_pub.publish(parameter_location_msg);
        }
        else if (!init_state.response.environment_surface_set)
        {
            ROS_INFO("Publishing environment surface\n");
            surface_pub.publish(surface_msg);
        }
        else if (!init_state.response.contact_sequence_set)
        {
            ROS_INFO("Publishing contact sequence\n");
            contact_sequence_pub.publish(contact_sequence_msg);
        }
        else if (!init_state.response.fully_initted)
        {
            ROS_INFO("Publishing init command\n");
            command_publisher.publish(galileo_cmd_msg);
            break;
        }

        ros::spinOnce();
        ros::Duration(0.3).sleep();
    }

    ros::spin();

    return 0;
}