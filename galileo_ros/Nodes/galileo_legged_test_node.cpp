#include "galileo_ros/GalileoLeggedRos.h"

std::vector<double> getX0(int nx, int q_index)
{
    galileo::opt::ConfigVector q0_vec = (galileo::opt::ConfigVector(19) << 0., 0., 0.339, 0., 0., 0., 1., 0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3).finished();

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

int main(int argc, char **argv)
{
    std::string solver_id = "go1_solver";
    std::string resources_location = argv[1];

    std::string model_location = resources_location + "/urdf/go1.urdf";
    std::string parameter_location = resources_location + "/SolverParameters/go1_solver_parameters.txt";

    ros::init(argc, argv);

    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    galileo::legged::GalileoLeggedRos go1_solver(nh, solver_id);

    // We will make calls from the same node, and see if they are read.
    galileo_ros::msg::ModelLocation model_location_msg;
    model_location_msg.model_location = model_location;
    model_location_msg.end_effector_names = {"FL_foot", "RL_foot", "FR_foot", "RR_foot"};
    galileo_ros::msg::ParameterLocation parameter_location_msg;
    parameter_location_msg.parameter_location = parameter_location;

    std::vector<int> knot_num = {180};
    std::vector<double> knot_time = {1.0};
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = {
        {0, 0, 0, 0}};

    galileo_ros::msg::ContactSequence contact_sequence_msg;
    contact_sequence_msg.knot_nums = knot_num;
    contact_sequence_msg.knot_times = knot_time;
    contact_sequence_msg.contact_surface_ids = contact_surfaces;

    ros::Publisher model_location_pub = nh->advertise<galileo_ros::msg::ModelLocation>(solver_id + "_model_location", 1);
    ros::Publisher parameter_location_pub = nh->advertise<galileo_ros::msg::ParameterLocation>(solver_id + "_parameter_location", 1);
    ros::Publisher contact_sequence_pub = nh->advertise<galileo_ros::msg::ContactSequence>(solver_id + "_contact_sequence", 1);

    model_location_pub.publish(model_location_msg);

    ros::spinOnce();

    parameter_location_pub.publish(parameter_location_msg);

    ros::spinOnce();

    contact_sequence_pub.publish(contact_sequence_msg);

    ros::spinOnce();

    ros::ServiceClient can_init_client = nh->serviceClient<std_srvs::Trigger>(solver_id + "can_init_service");
    std_srvs::Trigger can_init_srv;

    do
    {
        can_init_client.call(can_init_srv);
        ros::spinOnce();
    } while (!can_init_srv.response.success);

    galileo_ros::msg::GalileoCommand init_msg;

    std::vector<double> X0 = getX0(go1_solver.states()->nx, go1_solver.states()->q_index);

    init_msg.X_initial = X0;
    init_msg.X_goal = X0;

    ros::Publisher command_pub = nh->advertise<galileo_ros::msg::GalileoCommand>(solver_id + "_command", 1);
    init_pub.publish(init_msg);

    ros::spin();

    return 0;
}