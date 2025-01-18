#ifndef __galileo_estimator_measurement_model_base_hpp__
#define __galileo_estimator_measurement_model_base_hpp__

#include "galileo/estimator/fwd.hpp"

namespace galileo
{
    namespace estimator
    {

        template <typename _Derived>
        struct KalmanFilterBase;

        /**
         * @brief Abstract base class of all measurement models
         *
         * @param StateType The vector-type of the system state (usually some type derived from kalman::Vector)
         * @param MeasurementType The vector-type of the measurement (usually some type derived from kalman::Vector)
         * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
         */
        template <typename _Derived>
        class MeasurementModelBase : CRTP<_Derived>, public CovarianceBase<typename traits<_Derived>::Measurement>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            //! System state type
            using State = typename traits<_Derived>::State;

            //! System control input type
            using Measurement = typename traits<_Derived>::Measurement;

            /**
             * @brief Measurement model function h
             *
             * propagates the estimated measurement value
             * given the current state estimate x.
             */
            template <typename... Args>
            Measurement operator()(const State &x, Args &&...args) const
            {
                return derived().run(x, std::forward<Args>(args)...);
            }

        protected:
            template <typename>
            friend struct KalmanFilterBase;

            using CRTP<_Derived>::derived;

            GALILEO_DEFAULT_CONSTRUCTOR(MeasurementModelBase);
        };

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_measurement_model_base_hpp__