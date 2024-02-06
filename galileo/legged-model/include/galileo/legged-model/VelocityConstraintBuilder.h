#pragma once

#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {

            struct FootstepDefinition
            {
                double h_start;
                double h_end;
                double h_max;
            };

            casadi::SX createFootstepHeightFunction(casadi::SX &t, const FootstepDefinition &FS_def)
            {
                casadi::SX x0(2);
                x0(0) = if_else(t < 0.5, casadi::SX(FS_def.h_start), casadi::SX(FS_def.h_max));
                x0(1) = 0;
                casadi::SX xf(2);
                xf(0) = if_else(t >= 0.5, FS_def.h_end, FS_def.h_max);
                xf(1) = 0;

                casadi::SX H0(2, 2);
                H0(0, 0) = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
                H0(0, 1) = pow(t, 3) - 2 * pow(t, 2) + t;
                H0(1, 0) = 6 * pow(t, 2) - 6 * t;
                H0(1, 1) = 3 * pow(t, 2) - 4 * t + 1;

                casadi::SX Hf(2, 2);
                Hf(0, 0) = -2 * pow(t, 3) + 3 * pow(t, 2);
                Hf(0, 1) = pow(t, 3) - pow(t, 2);
                Hf(1, 0) = -6 * pow(t, 2) + 6 * t;
                Hf(1, 1) = 3 * pow(t, 2) - 2 * t;

                return casadi::SX::mtimes(H0, x0) + casadi::SX::mtimes(Hf, xf);
            }

            struct VelocityConstraintProblemData
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
                int num_knots;

                double max_footstep_offset_height;
                double corrector_kp;
                double following_leeway = 0;
            };

            template <class ProblemData>
            class VelocityConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                VelocityConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

                void BuildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data);

                void CreateBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
                {
                }

                void CreateFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
                {
                }
            };

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::BuildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data)
            {
                casadi::SXVector G_vec;
                casadi::SXVector upper_bound_vec;
                casadi::SXVector lower_bound_vec;

                auto t = problem_data.velocity_constraint_problem_data.t;

                contact::ContactSequence::CONTACT_SEQUENCE_ERROR error;

                contact::ContactMode mode = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(phase_index).mode;

                for (auto ee : problem_data.velocity_constraint_problem_data.robot_end_effectors)
                {
                    if (!mode[(*ee.second)])
                    {
                        std::cout << "EE " << ee.first << " is not in contact" << std::endl;
                        double liftoff_time = 0;
                        double touchdown_time = problem_data.velocity_constraint_problem_data.contact_sequence->dt();
                        FootstepDefinition footstep_definition;
                        footstep_definition.h_start = 0;
                        footstep_definition.h_end = 0;

                        // When did the EE break contact?
                        for (int i = phase_index; i >= 0; i--)
                        {
                            contact::ContactMode mode_i = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode;
                            // If the ee is in contact, we have found the phase where liftoff occurred.
                            if (mode_i.at(*ee.second))
                            {
                                int liftoff_index = i + 1;
                                problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(liftoff_index, liftoff_time, error);
                                auto surfaceID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee.second);
                                footstep_definition.h_start = (*problem_data.velocity_constraint_problem_data.environment_surfaces)[surfaceID].origin_z_offset;
                                break;
                            }
                        }

                        for (int i = phase_index; i < problem_data.velocity_constraint_problem_data.contact_sequence->num_phases() - 1; i++)
                        {
                            // If the ee is in contact, we have found the phase where touchdown occurred.
                            if (problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.at(*ee.second))
                            {
                                int touchdown_index = i;
                                problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(touchdown_index, touchdown_time, error);
                                auto surfaceID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee.second);
                                footstep_definition.h_end = (*problem_data.velocity_constraint_problem_data.environment_surfaces)[surfaceID].origin_z_offset;
                                break;
                                // if touchdown is not found, we will assume that the ee is in contact at the end of the horizon.
                                //  This is an odd and limiting assumption. Better behavior should be created.
                            }
                        }
                        footstep_definition.h_max = std::max(footstep_definition.h_start, footstep_definition.h_end) + problem_data.velocity_constraint_problem_data.max_footstep_offset_height;

                        // Add a baumgarte corrector height_velocity = (kp * (position(x) - desired_position) + desired_velocity)

                        // velocity_z - kp * (position_z(state)) = desired_velocity_z - kp * desired_position_z

                        // velocity_z - kp * (position_z(state)) < desired_velocity_z - kp * desired_position_z + leeway
                        // velocity_z - kp * (position_z(state)) > desired_velocity_z - kp * desired_position_z - leeway

                        casadi::Function t_in_range = casadi::Function("t_in_range", {t}, {(t - liftoff_time) / (touchdown_time - liftoff_time)});

                        casadi::SX footstep_height_expression = createFootstepHeightFunction(t_in_range(problem_data.velocity_constraint_problem_data.t).at(0), footstep_definition);
                        casadi::SX desired_height = footstep_height_expression(0);
                        casadi::SX desired_velocity = footstep_height_expression(1);

                        casadi::SX vel_bound = desired_velocity - problem_data.velocity_constraint_problem_data.corrector_kp * desired_height;

                        Eigen::Matrix<galileo::opt::ADScalar, 3, 1, 0> foot_pos = problem_data.velocity_constraint_problem_data.ad_data->oMf[ee.second->frame_id].translation();

                        casadi::SX cfoot_pos = casadi::SX(casadi::Sparsity::dense(foot_pos.rows(), 1));
                        pinocchio::casadi::copy(foot_pos, cfoot_pos);

                        auto foot_vel = pinocchio::getFrameVelocity(*(problem_data.velocity_constraint_problem_data.ad_model),
                                                                    *(problem_data.velocity_constraint_problem_data.ad_data),
                                                                    ee.second->frame_id,
                                                                    pinocchio::WORLD)
                                            .toVector();

                        casadi::SX cfoot_vel = casadi::SX(casadi::Sparsity::dense(foot_vel.rows(), 1));
                        pinocchio::casadi::copy(foot_vel, cfoot_vel);

                        casadi::SX vel_constraint_G = cfoot_vel - problem_data.velocity_constraint_problem_data.corrector_kp * cfoot_pos;

                        G_vec.push_back(vel_constraint_G);
                        lower_bound_vec.push_back(vel_bound - problem_data.velocity_constraint_problem_data.following_leeway);
                        upper_bound_vec.push_back(vel_bound + problem_data.velocity_constraint_problem_data.following_leeway);
                    }
                    else
                    {
                        // add velocity = 0 as a constraint
                        auto foot_vel = pinocchio::getFrameVelocity(*(problem_data.velocity_constraint_problem_data.ad_model),
                                                                    *(problem_data.velocity_constraint_problem_data.ad_data),
                                                                    ee.second->frame_id,
                                                                    pinocchio::WORLD)
                                            .toVector();

                        casadi::SX cfoot_vel = casadi::SX(casadi::Sparsity::dense(foot_vel.rows(), 1));
                        pinocchio::casadi::copy(foot_vel, cfoot_vel);

                        G_vec.push_back(cfoot_vel);
                        lower_bound_vec.push_back(casadi::SX::zeros(foot_vel.rows(), 1));
                        upper_bound_vec.push_back(casadi::SX::zeros(foot_vel.rows(), 1));
                    }
                }
                // replace u with an SX vector the size of (sum dof of ee in contact during this knot index)
                constraint_data.G = casadi::Function("G_Velocity", casadi::SXVector{problem_data.velocity_constraint_problem_data.x, problem_data.velocity_constraint_problem_data.u}, casadi::SXVector{casadi::SX::vertcat(G_vec)});
                std::cout << "G_velocity: " << constraint_data.G << std::endl;
                constraint_data.upper_bound = casadi::Function("upper_bound", casadi::SXVector{problem_data.velocity_constraint_problem_data.t}, casadi::SXVector{casadi::SX::vertcat(upper_bound_vec)});
                constraint_data.lower_bound = casadi::Function("lower_bound", casadi::SXVector{problem_data.velocity_constraint_problem_data.t}, casadi::SXVector{casadi::SX::vertcat(lower_bound_vec)});
            }
        }
    }
}