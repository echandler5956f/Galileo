#ifndef GALILEO_COMMON_UTILS_RESOURCE_FINDER_HPP
#define GALILEO_COMMON_UTILS_RESOURCE_FINDER_HPP

#include <string>
#include <filesystem>

#include "galileo/common/utils/resource_config.hpp"

namespace galileo
{
    namespace utils
    {
        /**
         * @brief Find the path to resources directory.
         *
         * This function tries to locate the resources directory in the following order:
         * 1. Source directory (for development builds)
         * 2. Installed location (for distribution builds)
         *
         * @return std::string Path to the resources directory
         * @throws std::runtime_error if resources directory cannot be found
         */
        inline std::string find_resources_dir()
        {
            // Try source directory first (for development)
            if (std::filesystem::exists(config::RESOURCES_DIR))
            {
                return config::RESOURCES_DIR;
            }

            // Try installed location (for distribution)
            if (std::filesystem::exists(config::INSTALLED_RESOURCES_DIR))
            {
                return config::INSTALLED_RESOURCES_DIR;
            }

            throw std::runtime_error("Cannot find Galileo resources directory. "
                                     "Expected locations: " +
                                     config::RESOURCES_DIR + " or " + config::INSTALLED_RESOURCES_DIR);
        }

        /**
         * @brief Get the path to a specific robot's URDF file.
         *
         * @param robot_name Name of the robot (e.g., "atlas", "go1", "huron")
         * @return std::string Full path to the robot's URDF file
         */
        inline std::string get_robot_urdf_path(const std::string &robot_name)
        {
            std::string resources_dir = find_resources_dir();
            return resources_dir + "/" + robot_name + "/urdf/" + robot_name + ".urdf";
        }

        /**
         * @brief Get the path to a specific robot's directory.
         *
         * @param robot_name Name of the robot (e.g., "atlas", "go1", "huron")
         * @return std::string Full path to the robot's directory
         */
        inline std::string get_robot_dir_path(const std::string &robot_name)
        {
            std::string resources_dir = find_resources_dir();
            return resources_dir + "/" + robot_name;
        }

        /**
         * @brief Get the path to a specific resource file.
         *
         * @param relative_path Relative path to the resource from the resources directory
         * @return std::string Full path to the resource file
         */
        inline std::string get_resource_path(const std::string &relative_path)
        {
            std::string resources_dir = find_resources_dir();
            return resources_dir + "/" + relative_path;
        }
    } // namespace utils
} // namespace galileo

#endif // GALILEO_COMMON_UTILS_RESOURCE_FINDER_HPP
