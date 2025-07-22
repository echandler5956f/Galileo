#ifndef __galileo_common_math_matrix_decomposition_hpp__
#define __galileo_common_math_matrix_decomposition_hpp__

#include "galileo/common/fwd.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <algorithm>

namespace galileo
{

    template <typename MatrixLike,
              bool value =
                  (Eigen::NumTraits<typename MatrixLike::Scalar>::IsInteger == 0)>
    struct pseudoInverseAlgo
    {
        using Scalar = typename MatrixLike::Scalar;
        using RealScalar = typename MatrixLike::RealScalar;

        static MatrixLike run(const Eigen::MatrixBase<MatrixLike> &a,
                              const RealScalar &epsilon)
        {
            using std::max;
            Eigen::JacobiSVD<MatrixLike> svd(a,
                                             Eigen::ComputeThinU | Eigen::ComputeThinV);
            RealScalar tolerance = epsilon *
                                   static_cast<Scalar>(max(a.cols(), a.rows())) *
                                   svd.singularValues().array().abs()(0);
            // FIX: Replace select() with a lambda function
            Eigen::Matrix<typename MatrixLike::Scalar, Eigen::Dynamic, 1>
                invSingularValues =
                    svd.singularValues().unaryExpr([&](const Scalar &x)
                                                   { return (x > tolerance) ? Scalar(1) / x : Scalar(0); });
            return svd.matrixV() * invSingularValues.asDiagonal() *
                   svd.matrixU().adjoint();
        }
    };

    template <typename MatrixLike>
    struct pseudoInverseAlgo<MatrixLike, false>
    {
        using Scalar = typename MatrixLike::Scalar;
        using RealScalar = typename MatrixLike::RealScalar;

        static MatrixLike run(const Eigen::MatrixBase<MatrixLike> &a,
                              const RealScalar &)
        {
            return Eigen::MatrixBase<MatrixLike>::Zero(a.rows(), a.cols());
        }
    };

    template <typename MatrixLike>
    MatrixLike pseudoInverse(
        const Eigen::MatrixBase<MatrixLike> &a,
        const typename MatrixLike::RealScalar &epsilon =
            Eigen::NumTraits<typename MatrixLike::Scalar>::dummy_precision())
    {
        return pseudoInverseAlgo<MatrixLike>::run(a, epsilon);
    }

} // namespace galileo

#endif // __galileo_common_math_matrix_decomposition_hpp__
