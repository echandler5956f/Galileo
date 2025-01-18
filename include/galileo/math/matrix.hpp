#ifndef __galileo_math_matrix_hpp__
#define __galileo_math_matrix_hpp__

#include "galileo/math/fwd.hpp"
#include <Eigen/Dense>

namespace galileo
{
    namespace math
    {

        /**
         * @brief Check if the input (square) matrix is symmetric
         *
         * @tparam _EigenDerived The EigenDerived type of the matrix
         * @param M The matrix to test for symmetry
         * @param eps The test tolerance
         * @return true if the matrix is symmetric
         * @return false if the matrix is not symmetric
         */
        template <typename _EigenDerived>
        bool isSymmetric(
            const Eigen::MatrixBase<_EigenDerived> &M,
            const typename _EigenDerived::Scalar eps = 1e-8)
        {
            return M.isApprox(M.transpose(), eps);
        }

        /**
         * @brief Enforce a (square) matrix to be symmetric
         *
         * @tparam _EigenDerived The EigenDerived type of the matrix
         * @param M The matrix to force symmetry on
         * @param eps The symmetry test tolerance
         * @return true if enforcing symmetry is successful, false otherwise
         */
        template <typename _EigenDerived>
        bool enforceSymmetric(
            Eigen::MatrixBase<_EigenDerived> &M,
            const typename _EigenDerived::Scalar eps = 1e-8)
        {
            using Scalar = typename _EigenDerived::Scalar;
            M = Scalar(0.5) * (M + M.transpose());
            return isSymmetric(M, eps);
        }

        /**
         * @brief Check if the input matrix is positive definite
         *
         * @tparam _EigenDerived The EigenDerived type of the matrix
         * @param M The matrix to test for positive definite
         * @param eps The test tolerance
         * @return true if the matrix is positive definite
         * @return false if the matrix is not positive definite
         */
        template <typename _EigenDerived>
        bool isPositiveDefinite(
            const Eigen::MatrixBase<_EigenDerived> &M,
            const typename _EigenDerived::Scalar eps = 1e-8)
        {
            Eigen::SelfAdjointEigenSolver<_EigenDerived> eigensolver(M);
            GALILEO_ASSERT(eigensolver.info() == Eigen::Success);
            if (eigensolver.info() == Eigen::Success)
            {
                // All eigenvalues must be >= 0:
                return (eigensolver.eigenvalues().array() >= eps).all();
            }
            return false;
        }

        /**
         * @brief Enforce a matrix to be positive definite
         *
         * @tparam _EigenDerived The EigenDerived type of the matrix
         * @param M The matrix to force for positive definite
         * @param eps The test tolerance
         * @return true if enforcing positive definite is successful, false otherwise
         */
        template <typename _EigenDerived>
        bool enforcePositiveDefinite(
            Eigen::MatrixBase<_EigenDerived> &M,
            const typename _EigenDerived::Scalar eps = 1e-8)
        {
            Eigen::SelfAdjointEigenSolver<_EigenDerived> eigensolver(M);
            GALILEO_ASSERT(eigensolver.info() == Eigen::Success);

            if (eigensolver.info() == Eigen::Success)
            {
                // All eigenvalues must be >= 0:
                using Scalar = typename _EigenDerived::Scalar;
                Scalar epsilon = eps;
                while ((eigensolver.eigenvalues().array() < eps).any())
                {
                    M.noalias() = eigensolver.eigenvectors() *
                                  eigensolver.eigenvalues().cwiseMax(epsilon).asDiagonal() *
                                  eigensolver.eigenvectors().transpose();
                    eigensolver.compute(M);
                    epsilon *= Scalar(10);
                }

                GALILEO_ASSERT(
                    isPositiveDefinite(M, eps),
                    "Failed to make matrix positive definite.");

                return epsilon != eps;
            }

            return false;
        }
    } // namespace math

    // template <typename Derived>
    // void test_eigen_matrix_base(Eigen::MatrixBase<Derived> &&s) {}

    // template <class, typename T>
    // struct is_eigen_matrix_impl : std::false_type
    // {
    // };

    // template <typename T>
    // struct is_eigen_matrix_impl<decltype(test_eigen_matrix_base(std::declval<T>())), T> : std::true_type
    // {
    // };

    // template <typename T>
    // struct is_eigen_matrix : is_eigen_matrix_impl<void, T>
    // {
    // };

    // template <typename T>
    // using enable_if_is_eigen_matrix = typename std::enable_if<is_eigen_matrix<T>::value>::type;

    // /**
    //  * @brief traits specialization for Eigen::Matrix
    //  */
    // template <typename Matrix>
    // struct traits<Matrix, enable_if_is_eigen_matrix<Matrix>>
    // {
    //     using Scalar = typename Matrix::Scalar;
    //     static constexpr auto Size = typename Matrix::RowsAtCompileTime;
    // };

} // namespace galileo

#endif // __galileo_math_matrix_hpp__