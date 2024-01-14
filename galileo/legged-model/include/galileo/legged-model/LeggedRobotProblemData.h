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
                casadi::SX t;
                casadi::SX x;
                casadi::SX u;
                std::shared_ptr<opt::LeggedRobotStates> states;
                std::shared_ptr<pinocchio::Model> model;
                contact::RobotEndEffectors robot_end_effectors;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;

                int num_knots;
            };
        }
    }
}