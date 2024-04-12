#include "galileo_ros/GalileoLeggedRos.h"

#include <galileo/tools/ReadFromFile.h>
#include <galileo/legged-model/LeggedModelHelpers.h>

void getProblemDataMessages(std::string urdf_name, std::string solver_parameter_file_name, std::string problem_parameter_file_name,
                            galileo_ros::RobotModel &robot_model_cmd,
                            galileo_ros::ParameterFileLocation &solver_parameter_location,
                            galileo_ros::ContactSequence &contact_sequence_cmd, galileo_ros::GalileoCommand &galileo_cmd_msg)
{

    std::vector<std::string> end_effectors;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;
    std::vector<double> q0;
    std::vector<double> qf;

    galileo::legged::helper::ReadProblemFromParameterFile(
        problem_parameter_file_name,
        end_effectors,
        knot_num,
        knot_time,
        contact_surfaces,
        q0,
        qf);

    std::string model_location = urdf_name;
    std::string solver_parameter_location_filename = solver_parameter_file_name;

    robot_model_cmd.model_file_location = model_location;

    robot_model_cmd.end_effector_names = end_effectors;

    solver_parameter_location.parameter_file_location = solver_parameter_location_filename;

    for (int phase_number = 0; phase_number < knot_num.size(); phase_number++)
    {
        galileo_ros::ContactPhase phase;
        phase.knot_num = knot_num[phase_number];
        phase.knot_time = knot_time[phase_number];
        phase.mode.contact_surface_ids = contact_surfaces[phase_number];
        contact_sequence_cmd.phases.push_back(phase);
    }

    galileo::legged::LeggedBody robot(model_location, end_effectors);

    std::vector<double> X0 = galileo::legged::helper::getXfromq(robot.si->nx, robot.si->q_index, q0);
    std::vector<double> Xf = galileo::legged::helper::getXfromq(robot.si->nx, robot.si->q_index, qf);

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
    std::string solver_parameter_file_name;
    std::string problem_parameter_file_name;

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/urdf_filename", urdf_file_name) ||
        !nh->getParam("galileo_ros/solver_parameters_location", solver_parameter_file_name) ||
        !nh->getParam("galileo_ros/problem_parameters_location", problem_parameter_file_name))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    galileo_ros::RobotModel model_location_msg;
    galileo_ros::ParameterFileLocation parameter_location_msg;
    galileo_ros::ContactSequence contact_sequence_msg;
    galileo_ros::GalileoCommand galileo_cmd_msg;
    galileo_ros::EnvironmentSurface surface_msg;

    std::vector<double> q0;
    std::vector<double> qf;
    getProblemDataMessages(urdf_file_name, solver_parameter_file_name, problem_parameter_file_name,
                           model_location_msg,
                           parameter_location_msg,
                           contact_sequence_msg, galileo_cmd_msg);

    ros::Publisher model_location_pub = nh->advertise<galileo_ros::RobotModel>(solver_id + "_model_location", 1);
    ros::Publisher parameter_location_pub = nh->advertise<galileo_ros::ParameterFileLocation>(solver_id + "_parameter_location", 1);
    ros::Publisher contact_sequence_pub = nh->advertise<galileo_ros::ContactSequence>(solver_id + "_contact_sequence", 1);
    ros::Publisher command_publisher = nh->advertise<galileo_ros::GalileoCommand>(solver_id + "_command", 1);
    ros::Publisher surface_pub = nh->advertise<galileo_ros::EnvironmentSurface>(solver_id + "_add_environment_surface", 1);
    ros::ServiceClient init_state_client = nh->serviceClient<galileo_ros::InitState>(solver_id + "_init_state_service");

    galileo_ros::InitState init_state;
    init_state.request.call = true;

    ros::spinOnce();
    while (!init_state_client.call(init_state))
    {
        ros::spinOnce();
        ros::Duration(0.1).sleep();
    }

    while (true)
    {
        init_state_client.call(init_state);
        if (!init_state.response.model_set)
        {
            model_location_pub.publish(model_location_msg);
        }
        else if (!init_state.response.solver_parameters_set)
        {
            parameter_location_pub.publish(parameter_location_msg);
        }
        else if (!init_state.response.environment_surface_set)
        {
            surface_pub.publish(surface_msg);
        }
        else if (!init_state.response.contact_sequence_set)
        {
            contact_sequence_pub.publish(contact_sequence_msg);
        }
        else if (!init_state.response.fully_initted)
        {
            command_publisher.publish(galileo_cmd_msg);
            break;
        }

        ros::spinOnce();
        ros::Duration(0.1).sleep();
    }

    ros::spin();

    return 0;
}