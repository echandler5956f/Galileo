#pragma once

#include "galileo/opt/PhaseSequence.h"

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Dummy mode class for the BasicSequence class.
         * 
         */
        class BasicMode
        {
            public:
                /**
                 * @brief Construct a new Basic Mode object.
                 * 
                 */
                BasicMode() {}

                /**
                 * @brief Destroy the Basic Mode object.
                 * 
                 */
                ~BasicMode() {}
        };

        /**
         * @brief Sequence class for basic phase sequences. Useful for simple single-phase problems, enabling usage with the TrajectoryOpt interface.
         * 
         */
        class BasicSequence : public PhaseSequence<BasicMode>
        {
            public:
                /**
                 * @brief Construct a new Basic Sequence object.
                 * 
                 */
                BasicSequence() {}

                /**
                 * @brief Destroy the Basic Sequence object.
                 * 
                 */
                ~BasicSequence() {}

                /**
                 * @brief Add a phase to the sequence.
                 * 
                 * @param mode The mode to add.
                 * @param knot_points The number of knot points.
                 * @param dt The time step.
                 * @param phase_dynamics The phase dynamics function.
                 * @return int The index of the added phase.
                 */
                int addPhase(const BasicMode &mode, int knot_points, double dt, casadi::Function phase_dynamics)
                {
                    return commonAddPhase(mode, knot_points, dt, phase_dynamics);
                }
        };
    }
}