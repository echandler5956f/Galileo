#include "galileo_ros/GalileoLeggedRos.h"
#include <fstream>
#include <sstream>
#include <filesystem>

std::vector<geometry_msgs::Point> loadSurfacePoints(std::string file_path)
{
    std::vector<geometry_msgs::Point> points;
    std::ifstream file;

    file.open(file_path);
    if (!file.is_open())
    {
        ROS_ERROR("Failed to open file %s", file_path.c_str());
        return points;
    }

    std::string line;
    std::cout << "Reading file " << file_path << std::endl;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z))
        {
            ROS_ERROR("Failed to read line %s", line.c_str());
            continue;
        }
        std::cout << "Read point " << x << " " << y << " " << z << std::endl;
        geometry_msgs::Point point;
        point.x = x;
        point.y = y;
        point.z = z;
        points.push_back( point );
    }

    file.close();

    return points;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "galileo_ros_environment_surface_node");
    std::shared_ptr<ros::NodeHandle> nh = std::make_shared<ros::NodeHandle>();

    std::string solver_id;
    std::string environment_surfaces_folder;

    if (!nh->getParam("galileo_ros/solver_id", solver_id) ||
        !nh->getParam("galileo_ros/environment_surfaces_location", environment_surfaces_folder) ){
        ROS_ERROR("Failed to get parameters");
        return 1;
        }

    std::vector<galileo_ros::EnvironmentSurface> surface_msgs;

    // Load the environment surfaces
    std::filesystem::recursive_directory_iterator it(environment_surfaces_folder);
    std::vector<std::string> files;
    for (auto &entry : it)
    {
        if (entry.is_regular_file())
        {
            files.push_back(entry.path().filename());
            ROS_INFO("Found file %s", entry.path().filename().c_str());
        }
    }



    for (auto file : files)
    {
        std::string file_path = environment_surfaces_folder + file;
        std::vector<geometry_msgs::Point> points = loadSurfacePoints(file_path);

        galileo_ros::EnvironmentSurface surface_msg;
        surface_msg.vertices = points;

        surface_msgs.push_back(surface_msg);
    }

    std::string topic = solver_id + "_add_environment_surface";
    
    ros::Publisher surface_pub = nh->advertise<galileo_ros::EnvironmentSurface>(topic, 1);

    ros::Duration(5).sleep();

    for (auto surface_msg : surface_msgs)
    {
        surface_pub.publish(surface_msg);
        ros::spinOnce();
    }

    ros::shutdown();
    return 0;
}