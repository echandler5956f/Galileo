#pragma once

#include "galileo/opt/PhaseSequence.h"

namespace galileo
{
    namespace opt
    {
        class BasicMode
        {
            public:
                BasicMode() {}
                ~BasicMode() {}
        };

        class BasicSequence : public PhaseSequence<BasicMode>
        {
            public:
                BasicSequence() {}
                ~BasicSequence() {}

                int addPhase(const BasicMode &mode, int knot_points, double dt, casadi::Function phase_dynamics)
                {
                    return commonAddPhase(mode, knot_points, dt, phase_dynamics);
                }
        };
    }
}