#include "variables/Constraint.h"

namespace galileo
{
    namespace model
    {
        namespace legged
        {
            class LeggedProblemData
            {
            public:
                LeggedProblemData();
                ~LeggedProblemData();

                variables::GeneralProblemData general_problem_data;
                model::legged::FrictionConeProblemData friction_cone_problem_data;
                model::legged::ContactConstraintProblemData contact_constraint_problem_data;
                model::legged::VelocityConstraintProblemData velocity_constraint_problem_data;
            };
        }
    }
}
