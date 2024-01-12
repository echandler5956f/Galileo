#pragma once

#include "galileo/legged-model/EnvironmentSurfaces.h"
#include <pinocchio/multibody/fwd.hpp>
#include <string>
#include <memory>
#include <map>

namespace galileo
{
    namespace legged
    {
        namespace contact
        {
            /**
             * @brief A struct for holding the end effector data of the robot.
             *
             */
            struct EndEffector
            {
                /**
                 * @brief The name of the end effector frame in pinocchio. Used as the key pair in global end effector maps.
                 *
                 */
                std::string frame_name;

                /**
                 * @brief The id of the frame in pinocchio.
                 *
                 */
                pinocchio::FrameIndex frame_id;

                /**
                 * @brief The ID the end effector is reffered to by in the vector of ee's.
                 *
                 */
                int local_ee_idx;

                /**
                 * @brief Is 6-DOF or is 3-DOF.
                 *
                 */
                bool is_6d;
            };

            /**
             * @brief We use a ptr so that there is only ever one instance of the end effector.
             *
             */
            typedef std::map<std::string, std::shared_ptr<EndEffector>> RobotEndEffectors;

            /**
             * @brief Maps the end effector name to the frame name.
             *
             */
            typedef std::map<std::string, bool> ContactCombination;
        }
    }
}