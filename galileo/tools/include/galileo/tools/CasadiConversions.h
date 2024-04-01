#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace tools
    {
        /**
         * @brief Convert an Eigen matrix to a Casadi matrix.
         *
         * @tparam Scalar Casadi matrix type
         * @param eigen_matrix The Eigen matrix to convert
         * @param casadi_matrix The Casadi matrix to store the result
         */
        template <typename Scalar>
        void eigenToCasadi(const Eigen::MatrixXd &eigen_matrix, Scalar &casadi_matrix);

        /**
         * @brief Convert a Casadi matrix to an Eigen matrix.
         *
         * @tparam Scalar Casadi matrix type
         * @param casadi_matrix The Casadi matrix to convert
         * @param eigen_matrix The Eigen matrix to store the result
         */
        template <typename Scalar>
        void casadiToEigen(const Scalar &casadi_matrix, Eigen::MatrixXd &eigen_matrix);

        /**
         * @brief Convert a vector to a Casadi matrix.
         *
         * @tparam Scalar Casadi matrix type
         * @param vec The vector to convert
         * @param rows The number of rows in the Casadi matrix
         * @param cols The number of columns in the Casadi matrix
         * @param casadi_matrix The Casadi matrix to store the result
         */
        template <typename Scalar>
        void vectorToCasadi(const std::vector<double>& vec, int rows, int cols, Scalar& casadi_matrix);

        /**
         * @brief Convert a Casadi matrix to a vector.
         *
         * @tparam Scalar Casadi matrix type
         * @param casadi_matrix The Casadi matrix to convert
         * @param vec The vector to store the result
         */
        template <typename Scalar>
        void casadiToVector(const Scalar& casadi_matrix, std::vector<double>& vec);

    }
}