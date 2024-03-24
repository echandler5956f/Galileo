#include <pinocchio/fwd.hpp>
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <visualization_msgs/Marker.h>
#include <std_msgs/Float64.h>

#include <galileo/legged-model/LeggedBody.h>

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

    // Get parameters from the parameter server
    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/urdf_filename", urdf_file_name))
    {
        ROS_ERROR("Failed to get parameters");
        return 1;
    }

    std::string end_effectors[] = {"FL_foot", "FR_foot", "RL_foot", "RR_foot"};
    int num_end_effectors = 4;

    // Load the URDF model to a leggedmodel
    galileo::legged::LeggedBody robot(urdf_file_name, num_end_effectors, end_effectors);

    ros::ServiceClient solution_client = nh->serviceClient<galileo_ros::SolutionRequest>(solver_id + "_get_solution");
    ros::Subscriber vis_time_sub = nh->subscribe("/visualization_time", 1, &visualizationTimeCallback);

    std::map<int, ros::Publisher> end_effector_position_pubs;
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
    double dt = 0.0123;

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

            for (int i = 0; i < num_times; i++)
            {
                std::vector<double> Xt = getSolutionAtTimeIdx(msg, i);

                std::vector<double> q_t(Xt.begin() + robot.si->q_index, Xt.begin() + robot.si->q_index + robot.si->nq);
                std::vector<double> u_t = getSolutionForceAtTimeIdx(msg, i);

                double ee_dof = u_t.size() / num_end_effectors;

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
            }

            for (auto &end_effector : robot.getEndEffectors())
            {
                visualization_msgs::Marker position_marker;
                position_marker.header.frame_id = "world";
                position_marker.header.stamp = ros::Time::now();
                position_marker.ns = end_effector.first + "_positions";
                position_marker.id = 0;
                position_marker.type = visualization_msgs::Marker::LINE_STRIP;
                position_marker.action = visualization_msgs::Marker::ADD;

                position_marker.scale.x = 0.01;
                position_marker.color.g = 0.30;
                position_marker.color.r = 0.30;
                position_marker.color.a = 1.0;
                position_marker.points = points[end_effector.first];
                end_effector_position_pubs[end_effector.first].publish(position_marker);
            }

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

            Eigen::VectorXd u = Eigen::Map<Eigen::VectorXd>(u_t.data(), u_t.size());

            int ee_dof = u_t.size() / num_end_effectors;
            auto robot_end_effectors = robot.getEndEffectors();

            int ee_local_index = 0;

            for (auto &end_effector : robot.getEndEffectors())
            {
                visualization_msgs::Marker force_marker;
                force_marker.header.frame_id = "world";
                force_marker.header.stamp = ros::Time::now();
                force_marker.ns = end_effector.first + "_forces";
                force_marker.id = 0;
                force_marker.type = visualization_msgs::Marker::ARROW;
                force_marker.action = visualization_msgs::Marker::ADD;

                force_marker.scale.x = 0.01;
                force_marker.scale.y = 0.01;
                force_marker.scale.z = 0.01;
                force_marker.color.b = 0.30;
                force_marker.color.r = 0.30;
                force_marker.color.a = 1.0;

                double force_scaling = 0.01;

                Eigen::Vector3d force = Eigen::Map<Eigen::VectorXd>(u_t.data() + ee_dof * ee_local_index, ee_dof).head(3);

                force = force * force_scaling;

                geometry_msgs::Point point = points[end_effector.first][index_of_time_closest_to_t];

                force_marker.points.push_back(point);

                geometry_msgs::Point termination_point;

                termination_point.x = point.x + force[0];
                termination_point.y = point.y + force[1];
                termination_point.z = point.z + force[2];

                std::cout << "Force: " << force.transpose() << std::endl;

                force_marker.points.push_back(termination_point);

                end_effector_force_publisher[end_effector.first].publish(force_marker);
                ee_local_index += 1;
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