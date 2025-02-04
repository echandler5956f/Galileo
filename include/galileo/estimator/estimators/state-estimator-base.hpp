#ifndef __galileo_estimator_state_estimator_base_hpp__
#define __galileo_estimator_state_estimator_base_hpp__

#include "galileo/estimator/fwd.hpp"

namespace galileo
{
    namespace estimator
    {
        template <typename _Scalar>
        class StateEstimatorBase
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        protected:
            Eigen::MatrixXd P; // Covariance matrix of the state estimate
            Eigen::MatrixXd R; // Covariance matrix of the measurement noise
            Eigen::MatrixXd Q; // Covariance matrix of the process noise

            size_t n_x; // Number of states
            size_t n_z; // Number of measurements
            size_t n_u; // Number of control inputs
            
        }; // class StateEstimatorBase

    } // namespace estimator

} // namespace galileo

#endif // __galileo_estimator_state_estimator_base_hpp__