#pragma once

#include <fstream>
#include <Eigen/Dense>

namespace galileo
{
    namespace tools
    {
        class MeshcatInterface
        {
        public:
            MeshcatInterface(std::string dir_);

            void WriteTimes(Eigen::VectorXd times, std::string file_name);

            void WriteJointPositions(Eigen::MatrixXd joint_positions, std::string file_name);

            void WriteMetadata(std::string urdf_location, Eigen::VectorXd q0, std::string file_name);

        protected:
            std::string dir;
        };
    }
}