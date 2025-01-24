#ifndef __galileo_math_matrix_hpp__
#define __galileo_math_matrix_hpp__

#include "galileo/math/fwd.hpp"
#include <Eigen/Dense>

namespace galileo
{
    namespace math
    {

        template <typename M1, typename M2>
        struct VConMat
        {
        private:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            static constexpr int R1 = Eigen::MatrixBase<M1>::RowsAtCompileTime;
            static constexpr int C1 = Eigen::MatrixBase<M1>::ColsAtCompileTime;
            static constexpr int R2 = Eigen::MatrixBase<M2>::RowsAtCompileTime;
            static constexpr int C2 = Eigen::MatrixBase<M2>::ColsAtCompileTime;

            // The new row count is either (R1+R2) if both are fixed, else Dynamic.
            static constexpr int Rows =
                (R1 != Eigen::Dynamic && R2 != Eigen::Dynamic)
                    ? (R1 + R2)
                    : Eigen::Dynamic;

            // If both have known, identical columns => keep that at compile time
            static constexpr bool sameFixedCols = (C1 != Eigen::Dynamic && C2 != Eigen::Dynamic && C1 == C2);
            static constexpr int Cols = sameFixedCols ? C1 : Eigen::Dynamic;

        public:
            // The underlying scalar type (e.g. double)
            using Scalar = typename M1::Scalar;

            // The resulting type
            using type = typename Eigen::Matrix<Scalar, Rows, Cols>;
        };

        template <typename V1, typename V2>
        struct VConVec
        {
        private:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            static constexpr int R1 = Eigen::MatrixBase<V1>::RowsAtCompileTime;
            static constexpr int R2 = Eigen::MatrixBase<V2>::RowsAtCompileTime;

            // The new row count is either (R1+R2) if both are fixed, else Dynamic.
            static constexpr int Rows =
                (R1 != Eigen::Dynamic && R2 != Eigen::Dynamic)
                    ? (R1 + R2)
                    : Eigen::Dynamic;

        public:
            // The underlying scalar type (e.g. double)
            using Scalar = typename V1::Scalar;

            // The resulting type
            using type = typename Eigen::Matrix<Scalar, Rows, 1>;
        };

        template <typename M1, typename M2>
        static typename VConMat<std::decay_t<M1>, std::decay_t<M2>>::type vertcat(M1 &&m1, M2 &&m2)
        {
            // 'std::decay_t<M1>' strips references and cv-qualifiers,
            // so if M1 is exactly some Eigen::Matrix<double,R,C> or expression,
            // we unify that type with VConMat.

            using M1Plain = std::decay_t<M1>;
            using M2Plain = std::decay_t<M2>;
            using ReturnType = typename VConMat<M1Plain, M2Plain>::type;

            // (1) Fully fixed in rows & cols
            if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                          ReturnType::ColsAtCompileTime != Eigen::Dynamic)
            {
                ReturnType res; // e.g. Matrix<double, R1+R2, C>
                res << m1, m2;  // one pass filling top/bottom
                return res;
            }
            // (2) Fixed cols, dynamic rows
            else if constexpr (ReturnType::RowsAtCompileTime == Eigen::Dynamic &&
                               ReturnType::ColsAtCompileTime != Eigen::Dynamic)
            {
                const int totalRows = m1.rows() + m2.rows();
                ReturnType res(totalRows, ReturnType::ColsAtCompileTime);
                res << m1, m2;
                return res;
            }
            // (3) Fixed rows, dynamic cols (unusual for vertical stacking, but included)
            else if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                               ReturnType::ColsAtCompileTime == Eigen::Dynamic)
            {
                const int totalCols = m1.cols(); // must match m2.cols() at runtime
                ReturnType res(ReturnType::RowsAtCompileTime, totalCols);
                res << m1, m2;
                return res;
            }
            // (4) Fully dynamic
            else
            {
                const int totalRows = m1.rows() + m2.rows();
                const int totalCols = m1.cols(); // must match m2.cols() at runtime
                ReturnType res(totalRows, totalCols);
                res << m1, m2;
                return res;
            }
        }

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