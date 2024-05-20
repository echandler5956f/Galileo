#include <pinocchio/fwd.hpp>
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <visualization_msgs/Marker.h>
#include <std_msgs/Float64.h>

#include <galileo/legged-model/LeggedBody.h>
#include <galileo/tools/ReadFromFile.h>

#include "galileo_ros/SolutionRequest.h"

std::vector<double> getSolutionAtTimeIdx(const galileo_ros::SolutionRequest &msg, int t_idx = 0)
{
    return std::vector<double>(msg.response.X_t_wrapped.begin() + msg.response.nx * t_idx, msg.response.X_t_wrapped.begin() + msg.response.nx * (1 + t_idx));
}

std::vector<double> getSolutionForceAtTimeIdx(const galileo_ros::SolutionRequest &msg, int t_idx = 0)
{
    return std::vector<double>(msg.response.U_t_wrapped.begin() + msg.response.nu * t_idx, msg.response.U_t_wrapped.begin() + msg.response.nu * (1 + t_idx));
}

double global_vis_time = 0;

void visualizationTimeCallback(const std_msgs::Float64::ConstPtr &msg)
{
    global_vis_time = msg->data;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "galileo_ros_legged_verbose_rviz_node");

    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    std::string solver_id;
    std::string urdf_file_name;
    std::string problem_parameter_file_name;

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/urdf_filename", urdf_file_name) ||
        !nh->getParam("galileo_ros/problem_parameters_location", problem_parameter_file_name))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    std::map<std::string, std::string> problem_parameters = galileo::tools::readFromFile(problem_parameter_file_name);

    std::vector<std::string> end_effectors_names = galileo::tools::readAsVector(problem_parameters["end_effector_names"]);
    std::vector<std::string> end_effector_types_str = galileo::tools::readAsVector(problem_parameters["end_effector_types"]);

    std::vector<galileo::legged::contact::EE_Types> end_effector_types;
    std::transform(end_effector_types_str.begin(), end_effector_types_str.end(), std::back_inserter(end_effector_types), [](const std::string &str)
                    { return static_cast<galileo::legged::contact::EE_Types>(std::stoi(str)); });

    // Load the URDF model to a leggedmodel
    galileo::legged::LeggedBody robot(urdf_file_name, end_effectors_names, end_effector_types);

    ros::ServiceClient solution_client = nh->serviceClient<galileo_ros::SolutionRequest>(solver_id + "_get_solution");
    ros::Subscriber vis_time_sub = nh->subscribe("/visualization_time", 1, &visualizationTimeCallback);

    std::map<int, ros::Publisher> end_effector_position_pubs;
    ros::Publisher com_traj_publisher = nh->advertise<visualization_msgs::Marker>("/com_traj", 1);
    ros::Publisher com_projected_publisher = nh->advertise<visualization_msgs::Marker>("/com_projected", 1);

    std::map<int, ros::Publisher> end_effector_force_publisher;

    for (auto &end_effector : robot.getEndEffectors())
    {
        std::string channel_name = std::string("/") + std::to_string(end_effector.first) + std::string("_positions");
        end_effector_position_pubs[end_effector.first] = nh->advertise<visualization_msgs::Marker>(channel_name, 1);

        channel_name = std::string("/") + std::to_string(end_effector.first) + std::string("_forces");
        end_effector_force_publisher[end_effector.first] = nh->advertise<visualization_msgs::Marker>(channel_name, 1);
    }

    galileo_ros::SolutionRequest msg;
    msg.request.times = {0};

    double horizon = 0.6;
    double dt = 0.0023;

    while (ros::ok())
    {
        if (solution_client.call(msg))
        {

            ros::spinOnce();
            if (!msg.response.solution_exists)
            {
                continue;
            }

            // get the positions of each end effector at each time and send them in a message
            int num_times = msg.response.times_evaluated.size();
            std::map<int, std::vector<geometry_msgs::Point>> points;

            std::vector<geometry_msgs::Point> com_traj;

            for (int i = 0; i < num_times; i++)
            {
                std::vector<double> Xt = getSolutionAtTimeIdx(msg, i);

                std::vector<double> q_t(Xt.begin() + robot.si->q_index, Xt.begin() + robot.si->q_index + robot.si->nq);

                Eigen::VectorXd q = Eigen::Map<Eigen::VectorXd>(q_t.data(), q_t.size());

                q.segment(3, 4) = q.segment(3, 4).normalized();

                pinocchio::framesForwardKinematics(robot.model, robot.data, q);
                pinocchio::updateFramePlacements(robot.model, robot.data);
                int ee_local_index = 0;
                for (auto &end_effector : robot.getEndEffectors())
                {
                    // GET POSITION
                    auto pin_point_ = robot.data.oMf[end_effector.first];

                    Eigen::Vector3d point_position = pin_point_.translation();
                    geometry_msgs::Point point;

                    point.x = point_position[0];
                    point.y = point_position[1];
                    point.z = point_position[2];

                    points[end_effector.first].push_back(point);
                }

                // GET COM POSITION
                Eigen::Vector3d com_pos = pinocchio::centerOfMass(robot.model, robot.data, q);
                geometry_msgs::Point com_point;
                com_point.x = com_pos[0];
                com_point.y = com_pos[1];
                com_point.z = com_pos[2];
                com_traj.push_back(com_point);
            }

            int ee_idx = 0;

            std::vector<std::vector<double>> force_colors{
                {38, 115, 126},
                {87, 187, 91},
                {139, 38, 53},
                {16, 29, 66}};

            for (auto &end_effector : robot.getEndEffectors())
            {
                visualization_msgs::Marker position_marker;
                position_marker.header.frame_id = "world";
                position_marker.header.stamp = ros::Time::now();
                position_marker.ns = end_effector.first + "_positions";
                position_marker.id = 0;
                position_marker.type = visualization_msgs::Marker::LINE_STRIP;
                position_marker.action = visualization_msgs::Marker::ADD;

                position_marker.scale.x = 0.005;
                position_marker.scale.y = 0.005;
                position_marker.scale.z = 0.005;

                position_marker.color.r = force_colors[ee_idx][0] / 255.0;
                position_marker.color.g = force_colors[ee_idx][1] / 255.0;
                position_marker.color.b = force_colors[ee_idx][2] / 255.0;

                position_marker.color.a = 1.0;
                position_marker.points = points[end_effector.first];
                end_effector_position_pubs[end_effector.first].publish(position_marker);

                ee_idx++;
            }

            visualization_msgs::Marker com_traj_marker;
            com_traj_marker.header.frame_id = "world";
            com_traj_marker.header.stamp = ros::Time::now();
            com_traj_marker.ns = "com_traj";
            com_traj_marker.id = 0;
            com_traj_marker.type = visualization_msgs::Marker::LINE_STRIP;
            com_traj_marker.action = visualization_msgs::Marker::ADD;

            com_traj_marker.scale.x = 0.015;
            com_traj_marker.scale.y = 0.015;
            com_traj_marker.scale.z = 0.015;

            com_traj_marker.color.r = 1.00;
            com_traj_marker.color.g = 0.40;
            com_traj_marker.color.b = 0.40;

            com_traj_marker.color.a = 1.0;
            com_traj_marker.points = com_traj;
            com_traj_publisher.publish(com_traj_marker);

            // get the forces of each end effector _now_ and send them in a message
            int index_of_time_closest_to_t = 0;

            for (int i = 0; i < num_times; i++)
            {
                if (msg.response.times_evaluated[i] > global_vis_time)
                {
                    index_of_time_closest_to_t = i;
                    break;
                }
            }

            std::vector<double> u_t = getSolutionForceAtTimeIdx(msg, index_of_time_closest_to_t);
            std::vector<double> X_t = getSolutionAtTimeIdx(msg, index_of_time_closest_to_t);

            std::vector<double> q_t(X_t.begin() + robot.si->q_index, X_t.begin() + robot.si->q_index + robot.si->nq);

            Eigen::VectorXd q = Eigen::Map<Eigen::VectorXd>(q_t.data(), q_t.size());

            pinocchio::framesForwardKinematics(robot.model, robot.data, q);
            pinocchio::updateFramePlacements(robot.model, robot.data);

            Eigen::VectorXd com_pos = pinocchio::centerOfMass(robot.model, robot.data, q);

            geometry_msgs::Point com_proj_point;

            com_proj_point.x = com_pos[0];
            com_proj_point.y = com_pos[1];
            com_proj_point.z = 0;

            visualization_msgs::Marker com_proj_marker;
            com_proj_marker.header.frame_id = "world";
            com_proj_marker.header.stamp = ros::Time::now();
            com_proj_marker.ns = "com_projected";
            com_proj_marker.id = 0;
            com_proj_marker.type = visualization_msgs::Marker::SPHERE;
            com_proj_marker.action = visualization_msgs::Marker::ADD;

            com_proj_marker.scale.x = 0.02;
            com_proj_marker.scale.y = 0.02;
            com_proj_marker.scale.z = 0.02;

            com_proj_marker.color.r = 0.50;
            com_proj_marker.color.g = 0.50;
            com_proj_marker.color.b = 0.50;

            com_proj_marker.color.a = 1.0;

            com_proj_marker.points.push_back(com_proj_point);
            com_proj_marker.points.push_back(com_proj_point);

            com_projected_publisher.publish(com_proj_marker);

            Eigen::VectorXd u = Eigen::Map<Eigen::VectorXd>(u_t.data(), u_t.size());
            auto robot_end_effectors = robot.getEndEffectors();
            int ee_local_start_index = 0;

            ee_idx = 0;
            for (auto &end_effector : robot.getEndEffectors())
            {
                visualization_msgs::Marker force_marker;
                force_marker.header.frame_id = "world";
                force_marker.header.stamp = ros::Time::now();
                force_marker.ns = end_effector.first + "_forces";
                force_marker.id = 0;
                force_marker.type = visualization_msgs::Marker::ARROW;
                force_marker.action = visualization_msgs::Marker::ADD;

                force_marker.scale.x = 0.03;
                force_marker.scale.y = 0.03;
                force_marker.scale.z = 0.01;

                force_marker.color.r = force_colors[ee_idx][0] / 255.0;
                force_marker.color.g = force_colors[ee_idx][1] / 255.0;
                force_marker.color.b = force_colors[ee_idx][2] / 255.0;

                force_marker.color.a = 1.0;

                double force_scaling = 0.0007;

                int ee_dof;

                if (end_effector.second->ee_type == galileo::legged::contact::EE_Types::NON_PREHENSILE_6DOF || end_effector.second->ee_type == galileo::legged::contact::EE_Types::PREHENSILE_6DOF)
                {
                    ee_dof = 6;
                }
                else if (end_effector.second->ee_type == galileo::legged::contact::EE_Types::NON_PREHENSILE_3DOF || end_effector.second->ee_type == galileo::legged::contact::EE_Types::PREHENSILE_3DOF)
                {
                    ee_dof = 3;
                }
                else
                {
                    throw std::runtime_error("End effector type not recognized.");
                }

                Eigen::Vector3d force = Eigen::Map<Eigen::VectorXd>(u_t.data() + ee_local_start_index, ee_dof).head(3);

                ee_local_start_index += ee_dof;

                force = force * force_scaling;

                geometry_msgs::Point point = points[end_effector.first][index_of_time_closest_to_t];

                force_marker.points.push_back(point);

                geometry_msgs::Point termination_point;

                termination_point.x = point.x + force[0];
                termination_point.y = point.y + force[1];
                termination_point.z = point.z + force[2];

                force_marker.points.push_back(termination_point);

                end_effector_force_publisher[end_effector.first].publish(force_marker);
                ee_idx++;
            }

            horizon = msg.response.solution_horizon;
            int N = horizon / dt;

            std::vector<double> call_times;
            for (int i = 0; i < N; i++)
            {
                call_times.push_back(i * horizon / N);
            }

            msg.request.times = call_times;
        }
    }

    return 0;
}