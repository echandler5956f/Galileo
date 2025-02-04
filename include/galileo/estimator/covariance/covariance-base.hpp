#ifndef __galileo_estimator_covariance_covariance_base_hpp__
#define __galileo_estimator_covariance_covariance_base_hpp__

#include "galileo/estimator/covariance/fwd.hpp"

namespace galileo
{
    namespace estimator
    {

        /**
         * @brief Base class for objects with Covariance
         *
         * @tparam StateType The state type
         */
        template <typename StateType>
        class CovarianceBase
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            /**
             * @brief Get covariance
             */
            const Covariance<StateType> &getCovariance() const
            {
                return P;
            }

            /**
             * @brief Set the covariance
             * @tparam Derived The Eigen-derived type of the input covariance
             * @param covariance The input covariance
             * @return true if the covariance is successfully set, false otherwise
             */
            template <typename EigenDerived>
            bool setCovariance(const Eigen::MatrixBase<EigenDerived> &covariance)
            {
                static_assert(
                    Eigen::MatrixBase<EigenDerived>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<EigenDerived>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");

                // Copy the covariance into the internal storage
                P = covariance.template cast<typename Covariance<StateType>::Scalar>();
                return true;
            }

            /**
             * @brief Get covariance (as square root)
             */
            CovarianceSquareRoot<StateType> getCovarianceSquareRoot() const
            {
                return CovarianceSquareRoot<StateType>(P);
            }

            /**
             * @brief Set Covariance using Square Root
             *
             * @param [in] covariance_square_root Lower triangular Matrix
             * representing the covariance square root (i.e. P = LLˆT).
             */
            template <typename EigenDerived>
            bool setCovarianceSquareRoot(const Eigen::MatrixBase<EigenDerived> &covariance_square_root)
            {
                static_assert(
                    Eigen::MatrixBase<EigenDerived>::RowsAtCompileTime == Covariance<StateType>::RowsAtCompileTime &&
                        Eigen::MatrixBase<EigenDerived>::ColsAtCompileTime == Covariance<StateType>::ColsAtCompileTime,
                    "Covariance matrix dimensions must match.");

                CovarianceSquareRoot<StateType> S;
                S.setL(covariance_square_root);
                return setCovariance(S.reconstructedMatrix());
            }

        protected:
            GALILEO_DEFAULT_CONSTRUCTOR(CovarianceBase);
            CovarianceBase(const Eigen::Ref<const Covariance<StateType>> &covariance)
                : P(covariance) {}

            //! Covariance
            Covariance<StateType> P = Covariance<StateType>::Identity() * 1e3;

        }; // class CovarianceBase

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_covariance_base_hpp__