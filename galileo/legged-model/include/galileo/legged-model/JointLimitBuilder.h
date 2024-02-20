#pragma once

#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {

            /**
             * @brief A struct for holding the data required to build the JointLimit Constraint.
             *
             */
            struct JointLimitProblemData
            {
                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<opt::LeggedRobotStates> states;
                std::shared_ptr<opt::ADModel> ad_model;
                std::shared_ptr<opt::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;
            };

            /**
             * A builder for the JointLimit Constraint.
             *
             * @tparam ProblemData must contain an instance of "JointLimitProblemData" named "joint_limit_problem_data"
             */
            template <class ProblemData>
            class JointLimitConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                JointLimitConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

            private:
                /**
                 * @brief Generate bounds for a vector of points.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "JointLimitProblemData" NAMED "joint_limit_problem_data"
                 * @param phase_index the index of the hase
                 * @param upper_bound Lower bound of the constraint at each point
                 * @param lower_bound Upper bound of the constraint at each point
                 */
                void createBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const;

                /**
                 * @brief Generate a function to evaluate each point
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "JointLimitProblemData" NAMED "joint_limit_problem_data"
                 * @param phase_index the index of the phase
                 * @param G A function that evaluates the constraint at each point
                 */
                void createFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const;
            };

            template <class ProblemData>
            void JointLimitConstraintBuilder<ProblemData>::createBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
            {
                // casadi::SX lb = vertcat(casadi::SXVector{
                //     -casadi::pi, -0.01, -1.2217, -casadi::pi / 2, -0.5236, -0.3491, -casadi::pi, -0.5236, -1.2217, -casadi::pi / 2, -0.5236, -0.3491});

                // casadi::SX ub = vertcat(casadi::SXVector{
                //     casadi::pi, 0.5236, 1.2217, casadi::pi / 2, 0.1745, 0.3491, casadi::pi, 0.01, 1.2217, casadi::pi / 2, 0.1745, 0.3491});

                casadi::SX lb = vertcat(casadi::SXVector{
                    -casadi::pi, -0.25, -1.2217, -casadi::pi / 2, -0.5236, -0.3491, -casadi::pi, -0.5236, -1.2217, -casadi::pi / 2, -0.5236, -0.3491});

                casadi::SX ub = vertcat(casadi::SXVector{
                    casadi::pi, 0.5236, 1.2217, casadi::pi / 2, 0.1745, 0.3491, casadi::pi, 0.25, 1.2217, casadi::pi / 2, 0.1745, 0.3491});

                lower_bound = casadi::Function("lower_bound", casadi::SXVector{problem_data.joint_limit_problem_data.t}, casadi::SXVector{lb});
                upper_bound = casadi::Function("upper_bound", casadi::SXVector{problem_data.joint_limit_problem_data.t}, casadi::SXVector{ub});
            }

            template <class ProblemData>
            void JointLimitConstraintBuilder<ProblemData>::createFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
            {
                G = casadi::Function("G_JointLimit", casadi::SXVector{problem_data.joint_limit_problem_data.x, problem_data.joint_limit_problem_data.u},
                                     casadi::SXVector{problem_data.joint_limit_problem_data.states->get_qj(problem_data.joint_limit_problem_data.x)});
                std::cout << "G_JointLimit: " << G << std::endl;
            }
        }
    }
}