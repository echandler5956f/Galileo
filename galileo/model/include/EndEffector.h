#pragma once

#include "EnvironmentSurfaces.h"
#include <pinocchio/multibody/fwd.hpp>
#include <string>
#include <memory>
#include <map>


namespace acro
{
    namespace contact
    {
        struct EndEffector
        {
            // The name of the end effector frame in pinocchio. Used as the key pair in global end effector maps.
            std::string frame_name;
            // The id of the frame in pinocchio
            pinocchio::FrameIndex frame_id;
            // The ID the end effector is reffered to by in the vector of ee's.
            int local_ee_idx;
            // Is 6-DOF or is 3-DOF
            bool is_6d;
        };

        // We use a ptr so that there is only ever one instance of the end effector.
        typedef std::map<std::string, std::shared_ptr<EndEffector>> RobotEndEffectors;

        typedef std::map<std::string, bool> ContactCombination;
    }
}