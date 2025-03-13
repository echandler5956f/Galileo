#ifndef __galileo_estimator_covariance_covariance_square_root_base_hpp__
#define __galileo_estimator_covariance_covariance_square_root_base_hpp__

#include "galileo/estimator/covariance/fwd.hpp"

namespace galileo
{
    namespace estimator
    {
        /**
         * @brief Base class for objects with Covariance as a square root
         *
         * @tparam StateType The state type
         */
        template <typename StateType>
        class CovarianceSquareRootBase
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            /**
             * @brief Get the reconstructed covariance matrix
             */
            Covariance<StateType> getCovariance() const
            {
                return S.reconstructedMatrix();
            }

            /**
             * @brief Set covariance as a covariance matrix
             * @param [in] covariance The input covariance
             * @return true if the covariance is successfully decomposed, false otherwise
             */
            template <typename CovarianceMatrixType>
            bool setCovariance(const Eigen::MatrixBase<CovarianceMatrixType> &covariance)
            {
                static_assert(
                    Eigen::MatrixBase<CovarianceMatrixType>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<CovarianceMatrixType>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");

                S.compute(covariance);
                return S.info() == Eigen::Success;
            }

            /**
             * @brief Get covariance as square root
             */
            const CovarianceSquareRoot<StateType> &getCovarianceSquareRoot() const
            {
                return S;
            }

            /**
             * @brief Set covariance using square root
             *
             * @param covariance_square_root Lower triangular matrix
             * representing the covariance square root (i.e. P = LLˆT).
             */
            template <typename LowerTriangularMatrixType>
            bool setCovarianceSquareRoot(const Eigen::MatrixBase<LowerTriangularMatrixType> &covariance_square_root)
            {
                static_assert(
                    Eigen::MatrixBase<LowerTriangularMatrixType>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<LowerTriangularMatrixType>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");

                CovarianceSquareRoot<StateType> S;
                S.setL(covariance_square_root);

                return true;
            }

        protected:
            GALILEO_DEFAULT_CONSTRUCTOR(CovarianceSquareRootBase);

            //! Covariance square root
            CovarianceSquareRoot<StateType> S = CovarianceSquareRoot<StateType>::Identity();

        }; // class CovarianceSquareRootBase

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_covariance_covariance_square_root_base_hpp__