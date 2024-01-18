#include "galileo/legged-model/FrictionConeConstraintBuilder.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {
            struct LeggedRobotProblemData
            {
                std::shared_ptr<opt::GeneralProblemData> gp_data;
                FrictionConeProblemData friction_cone_problem_data;
            };
        }
    }
}