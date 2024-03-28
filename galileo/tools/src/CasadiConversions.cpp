#include "galileo/tools/CasadiConversions.h"

namespace galileo
{
    namespace tools
    {
        template <typename Scalar>
        void eigenToCasadi(const Eigen::MatrixXd &eigen_matrix, Scalar &casadi_matrix)
        {
            size_t rows = eigen_matrix.rows();
            size_t cols = eigen_matrix.cols();

            casadi_matrix = Scalar::zeros(rows, cols);

            std::memcpy(casadi_matrix.ptr(), eigen_matrix.data(), sizeof(double) * rows * cols);
        }

        template void eigenToCasadi<casadi::DM>(const Eigen::MatrixXd &eigen_matrix, casadi::DM &casadi_matrix);
        // template void eigenToCasadi<casadi::SX>(const Eigen::MatrixXd &eigen_matrix, casadi::SX &casadi_matrix);
        // template void eigenToCasadi<casadi::MX>(const Eigen::MatrixXd &eigen_matrix, casadi::MX &casadi_matrix);

        template <typename Scalar>
        void casadiToEigen(const Scalar &casadi_matrix, Eigen::MatrixXd &eigen_matrix)
        {
            size_t rows = casadi_matrix.size1();
            size_t cols = casadi_matrix.size2();

            eigen_matrix.resize(rows, cols);
            eigen_matrix.setZero(rows, cols);

            std::memcpy(eigen_matrix.data(), casadi_matrix.ptr(), sizeof(double) * rows * cols);
        }

        template void casadiToEigen<casadi::DM>(const casadi::DM &casadi_matrix, Eigen::MatrixXd &eigen_matrix);
        // template void casadiToEigen<casadi::SX>(const casadi::SX &casadi_matrix, Eigen::MatrixXd &eigen_matrix);
        // template void casadiToEigen<casadi::MX>(const casadi::MX &casadi_matrix, Eigen::MatrixXd &eigen_matrix);

        template <typename Scalar>
        void vectorToCasadi(const std::vector<double> &vec, int rows, int cols, Scalar &casadi_matrix)
        {
            Eigen::MatrixXd intermediary_matrix = Eigen::Map<const Eigen::MatrixXd>(vec.data(), rows, cols);
            eigenToCasadi(intermediary_matrix, casadi_matrix);
        }

        template void vectorToCasadi<casadi::DM>(const std::vector<double> &vec, int rows, int cols, casadi::DM &casadi_matrix);
        // template void vectorToCasadi<casadi::SX>(const std::vector<double> &vec, int rows, int cols, casadi::SX &casadi_matrix);
        // template void vectorToCasadi<casadi::MX>(const std::vector<double> &vec, int rows, int cols, casadi::MX &casadi_matrix);

        template <typename Scalar>
        void casadiToVector(const Scalar &casadi_matrix, std::vector<double> &vec)
        {
            Eigen::MatrixXd intermediary_matrix;
            casadiToEigen(casadi_matrix, intermediary_matrix);
            vec = std::vector<double>(intermediary_matrix.data(), intermediary_matrix.data() + intermediary_matrix.size());
        }

        template void casadiToVector<casadi::DM>(const casadi::DM &casadi_matrix, std::vector<double> &vec);
        // template void casadiToVector<casadi::SX>(const casadi::SX &casadi_matrix, std::vector<double> &vec);
        // template void casadiToVector<casadi::MX>(const casadi::MX &casadi_matrix, std::vector<double> &vec);
    }
}