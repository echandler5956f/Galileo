#include "galileo/legged-model/FrictionConeConstraintBuilder.h"
#include "galileo/legged-model/VelocityConstraintBuilder.h"
#include "galileo/legged-model/ContactConstraintBuilder.h"

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
                VelocityConstraintProblemData velocity_constraint_problem_data;
                ContactConstraintProblemData contact_constraint_problem_data;                
            };
        }
    }
}