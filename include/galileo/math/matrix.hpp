#ifndef __galileo_math_matrix_hpp__
#define __galileo_math_matrix_hpp__

#include "galileo/math/fwd.hpp"
#include <Eigen/Dense>

namespace galileo
{
    namespace math
    {

        /**
         * @brief An alias for square matrix
         *
         * @tparam Scalar The scalar type
         * @tparam N The square dim
         * @tparam Options See Eigen
         * @tparam MaxRows See Eigen
         * @tparam MaxCols See Eigen
         *
         * @see Eigen::Matrix
         */
        template <typename Scalar, int N, int Options = 0, int MaxRows = N, int MaxCols = N>
        using SquareMatrix = Eigen::Matrix<Scalar, N, N, Options, MaxRows, MaxCols>;

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

                return epsilon != eps;
            }

            return false;
        }

        /**
         * @brief Check if the input matrix is a symmetric 
         * positive definite matrix
         *
         * @tparam MatrixType The MatrixType type of the matrix
         * @param M The matrix to test for symmetry and positive definiteness
         * @param eps The test tolerance
         * @return true is the matrix is SPD, false otherwise
         *
         * @see isSymmetric
         * @see isPositiveDefinite
         */
        template <typename MatrixType>
        static bool isSPD(
            const Eigen::MatrixBase<MatrixType> &M,
            const typename MatrixType::Scalar eps = 1e-8)
        {
            return isSymmetric(M, eps) && isPositiveDefinite(M, eps);
        }

        /**
         * @brief Enforce a matrix to be a symmetric 
         * positive definite matrix
         *
         * @tparam MatrixType The MatrixType type of the matrix
         * @param M The matrix to force as an SPD matrix
         * @param eps The test tolerance
         * @return true if enforcing SPD is successful, false otherwise
         *
         * @see enforceSymmetric
         * @see enforcePositiveDefinite
         */
        template <typename MatrixType>
        static bool enforceSPD(
            Eigen::MatrixBase<MatrixType> &M,
            const typename MatrixType::Scalar eps = 1e-8)
        {
            return enforceSymmetric(M, eps) && enforcePositiveDefinite(M, eps);
        }

        // This helper decays an arbitrary Eigen expression or matrix (Derived)
        // into a concrete Matrix type MType.
        //
        // If Derived is exactly MType and an rvalue, we can move.
        // Otherwise, we construct (i.e. evaluate) a new MType.
        template <typename Derived, typename MType>
        EIGEN_STRONG_INLINE MType toConcreteMatrix(Derived &&x)
        {
            // Check if the decayed type of x is exactly MType:
            if constexpr (std::is_same_v<std::decay_t<Derived>, MType>)
            {
                // We can move-construct MType if x is an rvalue of MType,
                // or copy-construct if it's an lvalue. std::forward picks the right one.
                return std::forward<Derived>(x);
            }
            else
            {
                // Evaluate the expression or convert from a different type
                // (e.g. from expression to Matrix, from block to Matrix, etc.)
                return MType(std::forward<Derived>(x));
            }
        }

    } // namespace math

} // namespace galileo

#endif // __galileo_math_matrix_hpp__