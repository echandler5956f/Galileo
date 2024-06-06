#pragma once

#include <vector>
#include <fstream>
#include <iostream>

#include "galileo/math/OrientationDefinition.h"
#include "galileo/legged-model/EndEffector.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"
#include "galileo/tools/ReadFromFile.h"
namespace galileo
{
    namespace legged
    {
        namespace helper
        {
            // Function to get the state (X) as a vector of doubles from the configuration (q). It is assumed the states are all 0 otherwise.
            std::vector<double> getXfromq(int nx, int q_index, std::vector<double> q);

            /**
             * @brief Reads the problem parameters from a file.
             *
             * @param problem_parameter_file_name Name of the file containing the problem parameters.
             * @param end_effector_names Vectors of strings containing the names of the end effectors.
             * @param end_effector_types Vector of integers containing the type of each end effector correspoding to the ee_names vector.
             * @param knot_num Vector of integers containing the number of knots for each phase.
             * @param knot_time Vector of doubles containing the time for each phase.
             * @param contact_surfaces Vector of phases, each phase is a vector of SurfaceID's, one for each end effector.
             * @param q0_vec Vector of doubles containing the initial configuration.
             * @param qf_vec Vector of doubles containing the target configuration.
             * @param orientation_def The internal representation of the orientation.
             */
            void ReadProblemFromParameterFile(std::string problem_parameter_file_name,
                                              std::vector<std::string> &end_effector_names,
                                              std::vector<contact::EE_Types> &end_effector_types,
                                              std::vector<int> &knot_num,
                                              std::vector<double> &knot_time,
                                              std::vector<std::vector<galileo::legged::environment::SurfaceID>> &contact_surfaces,
                                              std::vector<double> &q0_vec,
                                              std::vector<double> &qf_vec,
                                              math::OrientationDefinition &orientation_def);

            /**
             * @brief Get a Vector of binary masks from the Contact Surfaces. A contact surface value of NO_SURFACE is considered as not in contact.
             */
            std::vector<uint> getMaskVectorFromContactSurfaces(std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces);

        } // namespace helper
    } // namespace legged
} // namespace galileo