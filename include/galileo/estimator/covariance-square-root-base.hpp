#ifndef __galileo_estimator_covariance_square_root_base_hpp__
#define __galileo_estimator_covariance_square_root_base_hpp__

#include "galileo/estimator/fwd.hpp"

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
            template <typename _EigenDerived>
            bool setCovariance(const Eigen::MatrixBase<_EigenDerived> &covariance)
            {
                static_assert(
                    Eigen::MatrixBase<_EigenDerived>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<_EigenDerived>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");
                GALILEO_ASSERT(isCovariance(covariance), "CovarianceBase: Not a covariance matrix!");

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
            template <typename _EigenDerived>
            bool setCovarianceSquareRoot(const Eigen::MatrixBase<_EigenDerived> &covariance_square_root)
            {
                static_assert(
                    Eigen::MatrixBase<_EigenDerived>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<_EigenDerived>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");

                CovarianceSquareRoot<StateType> S;
                S.setL(covariance_square_root);

                GALILEO_ASSERT(isCovariance(S.reconstructedMatrix()));

                return true;
            }

        protected:
            GALILEO_DEFAULT_CONSTRUCTOR(CovarianceSquareRootBase);

            //! Covariance square root
            CovarianceSquareRoot<StateType> S = CovarianceSquareRoot<StateType>::Identity();
        };

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_covariance_square_root_base_hpp__