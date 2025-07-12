#ifndef GALILEO_TESTING_UTILS_RESOURCE_FINDER_HPP
#define GALILEO_TESTING_UTILS_RESOURCE_FINDER_HPP

#include <string>
#include <filesystem>

#include "utils/test_config.hpp"

namespace galileo
{
    namespace testing
    {
        /**
         * @brief Find the path to test resources directory.
         *
         * This function tries to locate the test resources directory in the following order:
         * 1. Source directory (for development builds)
         * 2. Installed location (for distribution builds)
         *
         * @return std::string Path to the resources directory
         * @throws std::runtime_error if resources directory cannot be found
         */
        inline std::string find_resources_dir()
        {
            // Try source directory first (for development)
            if (std::filesystem::exists(config::TEST_RESOURCES_DIR)) {
                return config::TEST_RESOURCES_DIR;
            }

            // Try installed location (for distribution)
            if (std::filesystem::exists(config::INSTALLED_RESOURCES_DIR)) {
                return config::INSTALLED_RESOURCES_DIR;
            }

            throw std::runtime_error("Cannot find Galileo test resources directory. "
                                   "Expected locations: " + config::TEST_RESOURCES_DIR +
                                   " or " + config::INSTALLED_RESOURCES_DIR);
        }

        /**
         * @brief Get the path to a specific robot's URDF file.
         *
         * @param robot_name Name of the robot (e.g., "atlas", "go1", "huron")
         * @return std::string Full path to the robot's URDF file
         */
        inline std::string get_robot_urdf_path(const std::string& robot_name)
        {
            std::string resources_dir = find_resources_dir();
            return resources_dir + "/" + robot_name + "/urdf/" + robot_name + ".urdf";
        }
    }
}

#endif // GALILEO_TESTING_UTILS_RESOURCE_FINDER_HPP
