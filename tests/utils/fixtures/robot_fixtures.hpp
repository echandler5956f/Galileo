#ifndef GALILEO_TESTING_UTILS_FIXTURES_ROBOT_FIXTURES_HPP
#define GALILEO_TESTING_UTILS_FIXTURES_ROBOT_FIXTURES_HPP

#include <string>
#include <vector>
#include <iostream>

#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/kinematics.hpp>

#include "galileo/multibody/robot-spec.hpp"
#include "galileo/multibody/states/implementations/multibody.hpp"
#include "galileo/multibody/actuations/implementations/floating-base.hpp"
#include "utils/resource_finder.hpp"

namespace galileo
{
    namespace testing
    {
        // A single, unified data structure to define a test robot case.
        // This is the core of the new scalable, data-driven fixture system.
        struct TestRobotData
        {
            std::string name;
            std::string urdf_path;
            bool is_floating_base;
            int nqj;
            int nvj;
        };

        // A single, central registry for all URDF-based test robots.
        // Adding a new robot to the test suite is as simple as adding a new entry here.
        inline std::vector<TestRobotData> get_test_robots()
        {
            return {
                {"Atlas", get_robot_urdf_path("atlas"), true, 30, 30},
                {"Go1", get_robot_urdf_path("go1"), true, 12, 12},
                {"Huron", get_robot_urdf_path("huron"), true, 12, 12},
            };
        }

        // A concrete RobotSpec that uses fully dynamic dimensions.
        // This allows us to use a single spec type for all URDF-based models.
        using DynamicSpec = galileo::RobotSpecTpl<
            galileo::BasicSpecTpl<double, double, 0>,
            detail::Dynamic, detail::Dynamic,
            detail::Dynamic, detail::Dynamic,
            detail::Dynamic,
            galileo::StateMultibodyTpl,
            galileo::ActuationFloatingBaseTpl>;

        // A fully-realized test fixture that holds the Pinocchio model/data
        // and the Galileo state/actuation models for a given robot.
        struct RobotModelFixture
        {
            RobotModelFixture(const TestRobotData &robot_data)
            {
                if (robot_data.is_floating_base)
                {
                    pinocchio::urdf::buildModel(robot_data.urdf_path, pinocchio::JointModelFreeFlyer(), model);
                }
                else
                {
                    pinocchio::urdf::buildModel(robot_data.urdf_path, model);
                }
                data = pinocchio::Data(model);

                state = std::make_shared<typename DynamicSpec::State_t>(spec, &model);
                actuation = std::make_shared<typename DynamicSpec::ActuationModel_t>(state);
            }

            DynamicSpec spec;
            pinocchio::Model model;
            pinocchio::Data data;

            std::shared_ptr<typename DynamicSpec::State_t> state;
            std::shared_ptr<typename DynamicSpec::ActuationModel_t> actuation;
        };
    }
}

#endif // GALILEO_TESTING_UTILS_FIXTURES_ROBOT_FIXTURES_HPP
