#ifndef __galileo_estimator_covariance_fwd_hpp__
#define __galileo_estimator_covariance_fwd_hpp__

#include "galileo/estimator/fwd.hpp"

namespace galileo
{
    namespace estimator
    {

        /**
         * @brief An alias for covariance matrix
         *
         * @tparam VectorType The vector type for which to generate a covariance
         * (usually a state or measurement type)
         *
         * @see SquareMatrix
         */
        template <typename VectorType>
        using Covariance = SquareMatrix<typename traits<VectorType>::Scalar, traits<VectorType>::Size>;

        /**
         * @brief An alias for covariance square root matrix
         * @param VectorType The vector type for which to generate a covariance
         * (usually a state or measurement type)
         *
         * @see Cholesky
         */
        template <typename VectorType>
        using CovarianceSquareRoot = math::Cholesky<Covariance<VectorType>>;

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_covariance_fwd_hpp__