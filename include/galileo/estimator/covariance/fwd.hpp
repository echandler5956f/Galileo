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

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_covariance_fwd_hpp__