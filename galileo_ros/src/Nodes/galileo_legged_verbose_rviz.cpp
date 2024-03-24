#include <pinocchio/fwd.hpp>
#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <visualization_msgs/Marker.h>

#include <galileo/legged-model/LeggedBody.h>

#include "galileo_ros/SolutionRequest.h"

std::vector<double> getSolutionAtTimeIdx(const galileo_ros::SolutionRequest &msg, int t_idx = 0)
{
    return std::vector<double>(msg.response.X_t_wrapped.begin() + msg.response.nx * t_idx, msg.response.X_t_wrapped.begin() + msg.response.nx * (1 + t_idx));
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
    // Load the URDF model to a leggedmodel
    galileo::legged::LeggedBody robot(urdf_file_name, 4, end_effectors);

    ros::ServiceClient solution_client = nh->serviceClient<galileo_ros::SolutionRequest>(solver_id + "_get_solution");

    ros::Publisher marker_pub = nh->advertise<visualization_msgs::Marker>("/end_effector_positions", 1);

    galileo_ros::SolutionRequest msg;
    msg.request.times = {0};

    double horizon = 0.4;
    double dt = 0.01;

    while (ros::ok())
    {
        if (solution_client.call(msg))
        {

            ros::spinOnce();
            if (!msg.response.solution_exists)
            {
                ros::spinOnce();
                continue;
            }

            // get the positions of each end effector at each time and send them in a message
            int num_times = msg.response.times_evaluated.size();
            std::vector<geometry_msgs::Point> points;

            for (int i = 0; i < num_times; i++)
            {
                std::vector<double> Xt = getSolutionAtTimeIdx(msg, i);
                std::vector<double> q_t(Xt.begin() + robot.si->q_index, Xt.begin() + robot.si->q_index + robot.si->nq);

                Eigen::VectorXd q = Eigen::Map<Eigen::VectorXd>(q_t.data(), q_t.size());

                q.segment(3, 4) = q.segment(3, 4).normalized();

                std::cout << "time " << msg.response.times_evaluated[i] << "; q: " << q.transpose() << std::endl;

                for (auto &end_effector : robot.getEndEffectors())
                {
                    pinocchio::forwardKinematics(robot.model, robot.data, q);
                    auto pin_point_ = robot.data.oMi[end_effector.first];

                    Eigen::Vector3d point_position = pin_point_.translation();
                    auto rot = pin_point_.rotation();

                    std::cout << "; q: " << q.transpose() << std::endl
                              << "rot: " << rot << std::endl;

                    std::cout << "time " << msg.response.times_evaluated[i] << "; q: " << q.transpose() << std::endl;

                    std::cout << "end effector " << end_effector.first << " FK: " << std::setprecision(10) << point_position << std::endl;
                    geometry_msgs::Point point;

                    point.x = point_position[0];
                    point.y = point_position[1];
                    point.z = point_position[2];

                    std::cout << "end effector " << end_effector.first << " position: " << point.x << " " << point.y << " " << point.z << std::endl;

                    points.push_back(point);
                }
            }

            visualization_msgs::Marker marker;
            marker.header.frame_id = "world";
            marker.header.stamp = ros::Time::now();
            marker.ns = "end_effector_positions";
            marker.id = 0;
            marker.type = visualization_msgs::Marker::POINTS;
            marker.action = visualization_msgs::Marker::ADD;
            marker.scale.x = 0.05;
            marker.scale.y = 0.05;
            marker.color.r = 1.0;
            marker.color.a = 1.0;
            marker.points = points;
            marker_pub.publish(marker);

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