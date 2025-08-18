#ifndef __galileo_common_math_matrix_decomposition_hpp__
#define __galileo_common_math_matrix_decomposition_hpp__

#include "galileo/common/fwd.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <algorithm>

namespace galileo
{

    template <typename MatrixLike, bool value = (Eigen::NumTraits<typename MatrixLike::Scalar>::IsInteger == 0)>
    struct pseudoInverseAlgo
    {
        using Scalar = typename MatrixLike::Scalar;
        using RealScalar = typename MatrixLike::RealScalar;
        // Pseudoinverse transposes dimensions: m×n matrix becomes n×m pseudoinverse
        using PseudoInverseType =
            Eigen::Matrix<Scalar, MatrixLike::ColsAtCompileTime, MatrixLike::RowsAtCompileTime, MatrixLike::Options>;

        static PseudoInverseType run(const Eigen::MatrixBase<MatrixLike> &a, const RealScalar &epsilon)
        {
            using std::max;
            // Convert to fully dynamic matrix to avoid Eigen's internal fixed-size workspace issues
            // while preserving the thin U/V optimization when possible
            using DynamicMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
            DynamicMatrix a_dynamic = a;

            constexpr unsigned int computationFlags = (MatrixLike::ColsAtCompileTime == Eigen::Dynamic)
                ? (Eigen::ComputeThinU | Eigen::ComputeThinV)
                : (Eigen::ComputeFullU | Eigen::ComputeFullV);
            Eigen::JacobiSVD<DynamicMatrix> svd(a_dynamic, computationFlags);
            RealScalar tolerance =
                epsilon * static_cast<Scalar>(max(a.cols(), a.rows())) * svd.singularValues().array().abs()(0);

            Eigen::Matrix<typename MatrixLike::Scalar, Eigen::Dynamic, 1> invSingularValues =
                svd.singularValues().unaryExpr([&](const Scalar &x)
                                               { return (x > tolerance) ? Scalar(1) / x : Scalar(0); });
            return svd.matrixV() * invSingularValues.asDiagonal() * svd.matrixU().adjoint();
        }
    };

    template <typename MatrixLike>
    struct pseudoInverseAlgo<MatrixLike, false>
    {
        using Scalar = typename MatrixLike::Scalar;
        using RealScalar = typename MatrixLike::RealScalar;
        // Pseudoinverse transposes dimensions: m×n matrix becomes n×m pseudoinverse
        using PseudoInverseType =
            Eigen::Matrix<Scalar, MatrixLike::ColsAtCompileTime, MatrixLike::RowsAtCompileTime, MatrixLike::Options>;

        static PseudoInverseType run(const Eigen::MatrixBase<MatrixLike> &a, const RealScalar &)
        {
            return PseudoInverseType::Zero(a.cols(), a.rows());
        }
    };

    template <typename MatrixLike>
    auto pseudoInverse(const Eigen::MatrixBase<MatrixLike> &a,
                       const typename MatrixLike::RealScalar &epsilon =
                           Eigen::NumTraits<typename MatrixLike::Scalar>::dummy_precision())
    {
        return pseudoInverseAlgo<MatrixLike>::run(a, epsilon);
    }

} // namespace galileo

#endif // __galileo_common_math_matrix_decomposition_hpp__
