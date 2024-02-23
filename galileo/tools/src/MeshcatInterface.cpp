#include "galileo/tools/MeshcatInterface.h"

namespace galileo
{
    namespace tools
    {
        MeshcatInterface::MeshcatInterface(std::string dir_)
        {
            this->dir = dir_;
        }

        void MeshcatInterface::WriteTimes(Eigen::VectorXd times, std::string file_name)
        {
            std::ofstream times_file(this->dir + file_name);
            if (times_file.is_open())
            {
                times_file << times.transpose().format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
                times_file.close();
            }
        }

        void MeshcatInterface::WriteJointPositions(Eigen::MatrixXd joint_positions, std::string file_name)
        {
            std::ofstream joint_positions_file(this->dir + file_name);
            if (joint_positions_file.is_open())
            {
                joint_positions_file << joint_positions.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
                joint_positions_file.close();
            }
        }

        void MeshcatInterface::WriteMetadata(std::string urdf_location, Eigen::VectorXd q0, std::string file_name)
        {
            std::ofstream file(this->dir + file_name);
            if (file.is_open())
            {
                file << "urdf location: " << urdf_location << "\n";

                file << "q0: ";
                for (int i = 0; i < q0.size(); ++i)
                {
                    file << q0[i];
                    if (i != q0.size() - 1) // not the last element
                    {
                        file << ", ";
                    }
                }
                file << "\n";

                file.close();
            }
        }
    }
}