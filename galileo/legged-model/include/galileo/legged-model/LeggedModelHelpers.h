#pragma once

#include <vector>
#include <fstream>
#include <iostream>

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
             * @brief Reads a vector of "phases" comprised of contact surfaces.
             *
             * @param contact_surfaces_str Vector of strings, Element "i" corresponds to the "i"th phase and is formatted as "(surface_for_end_effector_1,surface_for_end_effector_2,...,surface_for_end_effector_n)"
             * @return std::vector<std::vector<galileo::legged::environment::SurfaceID>> Vector of phases, each phase is a vector of SurfaceID's, one for each end effector
             *
             */
            std::vector<std::vector<galileo::legged::environment::SurfaceID>> readContactSurfaces(std::vector<std::string> contact_surfaces_str);

            /**
             * @brief Reads the problem parameters from a file.
             *
             * @param problem_parameter_file_name Name of the file containing the problem parameters.
             * @param end_effector_names Vector of strings containing the names of the end effectors.
             * @param knot_num Vector of integers containing the number of knots for each phase.
             * @param knot_time Vector of doubles containing the time for each phase.
             * @param contact_surfaces Vector of phases, each phase is a vector of SurfaceID's, one for each end effector.
             *
             * @param q0_vec Vector of doubles containing the initial configuration.
             * @param qf_vec Vector of doubles containing the target configuration.
             */
            void ReadProblemFromParameterFile(std::string problem_parameter_file_name,
                                              std::vector<std::string> &end_effector_names,
                                              std::vector<int> &knot_num,
                                              std::vector<double> &knot_time,
                                              std::vector<std::vector<galileo::legged::environment::SurfaceID>> &contact_surfaces,
                                              std::vector<double> &q0_vec,
                                              std::vector<double> &qf_vec);

            /**
             * @brief Get a Vector of binary masks from the Contact Surfaces. A contact surface value of NO_SURFACE is considered as not in contact.
             */
            std::vector<uint> getMaskVectorFromContactSurfaces(std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces);

        } // namespace helper
    }     // namespace legged
} // namespace galileo