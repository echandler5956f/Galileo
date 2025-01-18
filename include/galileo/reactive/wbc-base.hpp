#ifndef __galileo_reactive_wbc_base_hpp__
#define __galileo_reactive_wbc_base_hpp__

#include "galileo/reactive/fwd.hpp"
#include "galileo/reactive/task.hpp"

namespace galileo
{
    namespace reactive
    {
        // Decision Variables: x = [\dot u^T, F^T, \tau^T]^T
        class WbcBase
        {
            using Vector6 = Eigen::Matrix<scalar_t, 6, 1>;
            using Matrix6 = Eigen::Matrix<scalar_t, 6, 6>;

        public:
            WbcBase(const PinocchioInterface &pinocchioInterface, CentroidalModelInfo info, const PinocchioEndEffectorKinematics &eeKinematics);

            virtual vector_t update(const vector_t &stateDesired, const vector_t &inputDesired, const vector_t &rbdStateMeasured, size_t mode,
                                    scalar_t period);

        protected:
            void updateMeasured(const vector_t &rbdStateMeasured);
            void updateDesired(const vector_t &stateDesired, const vector_t &inputDesired);

            size_t getNumDecisionVars() const { return numDecisionVars_; }

            size_t numDecisionVars_;
            PinocchioInterface pinocchioInterfaceMeasured_, pinocchioInterfaceDesired_;
            CentroidalModelInfo info_;

            std::unique_ptr<PinocchioEndEffectorKinematics> eeKinematics_;
            CentroidalModelPinocchioMapping mapping_;

            vector_t qMeasured_, vMeasured_, inputLast_;
            matrix_t j_, dj_;
            contact_flag_t contactFlag_{};
            size_t numContacts_{};

            // Task Parameters:
            vector_t torqueLimits_;
            scalar_t frictionCoeff_{}, swingKp_{}, swingKd_{};
        };

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_wbc_base_hpp__