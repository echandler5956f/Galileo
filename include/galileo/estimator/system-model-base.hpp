#ifndef __galileo_estimator_system_model_base_hpp__
#define __galileo_estimator_system_model_base_hpp__

#include "galileo/estimator/fwd.hpp"

namespace galileo
{
    namespace estimator
    {
        template <typename _Derived>
        struct KalmanFilterBase;

        /**
         * @brief Abstract base class of all system models
         *
         * @param StateType The vector-type of the system state (usually some type derived from kalman::Vector)
         * @param ControlType The vector-type of the control input (usually some type derived from kalman::Vector)
         * @param CovarianceBase The class template used for covariance storage (must be either StandardBase or SquareRootBase)
         */
        template <typename _Derived>
        class SystemModelBase : CRTP<_Derived>, public CovarianceBase<typename traits<_Derived>::Control>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            //! System state type
            using State = typename traits<_Derived>::State;

            //! State Scalar type
            using Scalar = typename traits<State>::Scalar;

            //! System control input type
            using Control = typename traits<_Derived>::Control;

            /**
             * @brief State transition function f
             *
             * Computes the propagateed system state in the next timestep given
             * the current state x and the control input u
             *
             * @return The propagated system state
             */
            template <typename... Args>
            State operator()(Args &&...args) const
            {
                return derived().run(std::forward<Args>(args)...);
            }

        protected:
            using Base = CovarianceBase<typename traits<_Derived>::Control>;

            template <typename>
            friend class KalmanFilterBase;

            using CRTP<_Derived>::derived;
            using Base::Base;

            GALILEO_DEFAULT_CONSTRUCTOR(SystemModelBase);
        };

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_system_model_base_hpp__