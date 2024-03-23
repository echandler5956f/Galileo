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

void getProblemDataMessages(std::string resources_location, std::string urdf_name, std::string parameter_file_name,
                            galileo_ros::RobotModel &robot_model_cmd,
                            galileo_ros::ParameterFileLocation &parameter_location_cmd,
                            galileo_ros::ContactSequence &contact_sequence_cmd, galileo_ros::GalileoCommand &galileo_cmd_msg)
{
    std::string model_location = resources_location + "/urdf/" + urdf_name;
    std::string parameter_location = resources_location + "/SolverParameters/" + parameter_file_name;

    robot_model_cmd.model_file_location = model_location;
    robot_model_cmd.end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};

    parameter_location_cmd.parameter_file_location = parameter_location;

    std::vector<int> knot_num = {180};
    std::vector<double> knot_time = {1.0};
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
        {0, 0, 0, 0}};

    for (int phase_number = 0; phase_number < knot_num.size(); phase_number++)
    {
        galileo_ros::ContactPhase phase;
        phase.knot_num = knot_num[phase_number];
        phase.knot_time = knot_time[phase_number];
        phase.contact_surface_ids = contact_surfaces[phase_number];
        contact_sequence_cmd.phases.push_back(phase);
    }

    std::string end_effectors[] = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};

    galileo::legged::LeggedBody robot(model_location, 4, end_effectors);

    std::vector<double> X0 = getX0(robot.si->nx, robot.si->q_index);

    galileo_cmd_msg.command_type = "init";

    galileo_cmd_msg.X_initial = X0;
    galileo_cmd_msg.X_goal = X0;
}

int main(int argc, char **argv)
{
    // ros::start();
    // ros::Rate loop_rate(10);

    if (argc != 5)
    {
        std::cerr << "Usage: rosrun galileo_ros galileo_legged_test_node <solver_id> <resources_location> <urdf_file_name> <parameter_file_name>" << std::endl;
        std::cerr << "<resources_location>/urdf/<urdf_file_name> and <resources_location>/SolverParameters/<parameter_file_name> must exist" << std::endl;
        return 1;
    }

    std::string solver_id = argv[1];
    std::string resources_location = argv[2];
    std::string urdf_file_name = argv[3];
    std::string parameter_file_name = argv[4];

    ros::init(argc, argv, solver_id + "_test_node");
    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    ROS_INFO("Creating the solver object");

    // galileo::legged::GalileoLeggedRos go1_solver(nh, solver_id);

    galileo_ros::RobotModel model_location_msg;
    galileo_ros::ParameterFileLocation parameter_location_msg;
    galileo_ros::ContactSequence contact_sequence_msg;
    galileo_ros::GalileoCommand galileo_cmd_msg;

    ROS_INFO("Generating problem data messages");

    getProblemDataMessages(resources_location, urdf_file_name, parameter_file_name,
                           model_location_msg,
                           parameter_location_msg,
                           contact_sequence_msg, galileo_cmd_msg);

    ROS_INFO("Publishing model location, parameter location and contact sequence");

    ros::Publisher model_location_pub = nh->advertise<galileo_ros::RobotModel>(solver_id + "_model_location", 1);
    ros::Publisher parameter_location_pub = nh->advertise<galileo_ros::ParameterFileLocation>(solver_id + "_parameter_location", 1);
    ros::Publisher contact_sequence_pub = nh->advertise<galileo_ros::ContactSequence>(solver_id + "_contact_sequence", 1);
    ros::Publisher command_publisher = nh->advertise<galileo_ros::GalileoCommand>(solver_id + "_command", 1);
    ros::ServiceClient init_state_client = nh->serviceClient<galileo_ros::InitState>(solver_id + "_init_state_service");

    galileo_ros::InitState init_state;
    init_state.request.call = true;

    ROS_INFO("calling init state");
    ros::spinOnce();
    while (!init_state_client.call(init_state))
    {
        ROS_INFO("still init state");
        ros::spinOnce();
        ros::Duration(0.3).sleep();
    }

    ROS_INFO("called init state");

    while (true)
    {
        init_state_client.call(init_state);
        ROS_INFO("model set: %d, solver parameters set: %d, contact sequence set: %d, fully initted: %d",
                 init_state.response.model_set,
                 init_state.response.solver_parameters_set,
                 init_state.response.contact_sequence_set,
                 init_state.response.fully_initted);

        if (!init_state.response.model_set)
        {

            ROS_INFO("Publishing model location");
            model_location_pub.publish(model_location_msg);
        }
        else if (!init_state.response.solver_parameters_set)
        {
            ROS_INFO("Publishing parameter location");
            parameter_location_pub.publish(parameter_location_msg);
        }
        else if (!init_state.response.contact_sequence_set)
        {
            ROS_INFO("Publishing contact sequence");
            contact_sequence_pub.publish(contact_sequence_msg);
        }
        else if (!init_state.response.fully_initted)
        {
            ROS_INFO("Publishing init command");
            command_publisher.publish(galileo_cmd_msg);
            break;
        }

        ros::spinOnce();
        ros::Duration(0.3).sleep();
    }

    ros::spin();

    return 0;
}