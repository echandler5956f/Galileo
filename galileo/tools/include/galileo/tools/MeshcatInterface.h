#pragma once

#include <fstream>
#include <Eigen/Dense>

namespace galileo
{
    namespace tools
    {
        /**
         * @brief Simple interface to write solutions to file, which can be read in Python to use Meshcat
         * 
         */
        class MeshcatInterface
        {
        public:
            /**
             * @brief Construct a new MeshcatInterface object.
             * 
             * @param dir_ Directory where the files will be saved.
             */
            MeshcatInterface(std::string dir_);

            /**
             * @brief Write the times to a file.
             * 
             * @param times The times to be written.
             * @param file_name The name of the file to be written.
             */
            void WriteTimes(Eigen::VectorXd times, std::string file_name);

            /**
             * @brief Write the joint positions to a file.
             * 
             * @param joint_positions The joint positions to be written.
             * @param file_name The name of the file to be written.
             */
            void WriteJointPositions(Eigen::MatrixXd joint_positions, std::string file_name);

            /**
             * @brief Write the urdf location to a file as metadata.
             * 
             * @param urdf_location The location of the urdf file.
             * @param file_name The name of the file to be written.
             */
            void WriteMetadata(std::string urdf_location, std::string file_name);

        protected:
            std::string dir;
        };
    }
}