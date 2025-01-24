#ifndef __galileo_estimator_types_hpp__
#define __galileo_estimator_types_hpp__

#include "galileo/estimator/fwd.hpp"

namespace galileo
{
    namespace estimator
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
         * @brief An alias for covariance matrix
         *
         * @tparam Type The vector type for which to generate a covariance
         * (usually a state or measurement type)
         *
         * @see SquareMatrix
         */
        template <typename Type>
        using Covariance = SquareMatrix<typename traits<Type>::Scalar, traits<Type>::Size>;

        /**
         * @brief An alias for covariance square root matrix
         * @param Type The vector type for which to generate a covariance
         * (usually a state or measurement type)
         *
         * @see Cholesky
         */
        template <typename Type>
        using CovarianceSquareRoot = math::Cholesky<Covariance<Type>>;

        /**
         * @brief An alias for the Kalman Gain matrix
         * @param State The system state type
         * @param Measurement The measurement type
         */
        template <typename State, typename Measurement, int Options = 0, int MaxRows = traits<State>::Size, int MaxCols = Measurement::RowsAtCompileTime>
        using KalmanGain = Eigen::Matrix<typename traits<State>::Scalar, traits<State>::Size, Measurement::RowsAtCompileTime, Options, MaxRows, MaxCols>;

        /**
         * @brief An alias for the Jacobian matrix J_A_B
         *
         * @tparam A The 'input' vector type (usually a state or measurement type)
         * @tparam B The 'output' vector type (usually a tangent or measurement type)
         * @tparam Options see Eigen
         * @tparam MaxRows see Eigen
         * @tparam MaxCols see Eigen
         *
         * @see Eigen::Matrix
         */
        template <typename A, typename B, int Options = 0, int MaxRows = traits<A>::Size, int MaxCols = traits<B>::Size>
        using Jacobian = Eigen::Matrix<typename traits<A>::Scalar, traits<A>::Size, traits<B>::Size, Options, MaxRows, MaxCols>;

        /**
         * @brief Check if the input matrix is a covariance matrix
         * (symmetric positive definite matrix)
         *
         * @tparam EigenDerived The EigenDerived type of the matrix
         * @param M The matrix to test for covariance
         * @param eps The test tolerance
         * @return true is the matrix is a covariance, false otherwise
         *
         * @see isSymmetric
         * @see isPositiveDefinite
         */
        template <typename EigenDerived>
        static bool isCovariance(
            const Eigen::MatrixBase<EigenDerived> &M,
            const typename EigenDerived::Scalar eps = 1e-8)
        {
            return math::isSymmetric(M, eps) && math::isPositiveDefinite(M, eps);
        }

        /**
         * @brief Enforce a matrix to be a covariance matrix
         * (symmetric positive definite matrix)
         *
         * @tparam EigenDerived  The EigenDerived type of the matrix
         * @param M The matrix to force as a covariance matrix
         * @param eps The test tolerance
         * @return true if enforcing covariance is successful, false otherwise
         *
         * @see enforceSymmetric
         * @see enforcePositiveDefinite
         */
        template <typename EigenDerived>
        static bool enforceCovariance(
            Eigen::MatrixBase<EigenDerived> &M,
            const typename EigenDerived::Scalar eps = 1e-8)
        {
            return math::enforceSymmetric(M, eps) && math::enforcePositiveDefinite(M, eps);
        }

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_types_hpp__