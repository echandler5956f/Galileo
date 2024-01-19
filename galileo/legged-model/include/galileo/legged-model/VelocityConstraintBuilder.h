#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"
#include "galileo/legged-model/LeggedRobotStates.h"
#include "galileo/legged-model/FootstepTrajectoryGenerator.h"

namespace
{
    casadi::Function createFootstepHeightFunction(casadi::SX &t, const FootstepDefinition &FS_def)
    {
        casadi::MX x0(2);
        x0(0) = casadi::ifelse(t < 0.5, FS_def.h_start, FS_def.h_max);
        x0(1) = 0;
        casadi::MX xf(2);
        xf(0) = casadi::ifelse(t >= 0.5, FS_def.h_end, FS_def.h_max);
        xf(1) = 0;

        casadi::MX H0(2, 2);
        H0(0, 0) = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
        H0(0, 1) = pow(t, 3) - 2 * pow(t, 2) + t;
        H0(1, 0) = 6 * pow(t, 2) - 6 * t;
        H0(1, 1) = 3 * pow(t, 2) - 4 * t + 1;

        casadi::MX Hf(2, 2);
        Hf(0, 0) = -2 * pow(t, 3) + 3 * pow(t, 2);
        Hf(0, 1) = pow(t, 3) - pow(t, 2);
        Hf(1, 0) = -6 * pow(t, 2) + 6 * t;
        Hf(1, 1) = 3 * pow(t, 2) - 2 * t;

        casadi::Function footstep_height_function = casadi::Function("footstep_z_function", {t}, {H0 * x0 + Hf * xf});

        return footstep_height_function;
    }

    struct FootstepDefinition
    {
        double h_start;
        double h_end;
        double h_max;
    };

}
namespace galileo
{
    namespace legged
    {
        namespace constraints
        {
            struct VelocityConstraintProblemData
            {

                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<opt::States> states;
                std::shared_ptr<opt::Model> model;
                std::shared_ptr<opt::Data> data;
                std::shared_ptr<opt::ADModel> ad_model;
                std::shared_ptr<opt::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x; // this needs to be initialized to casadi::SX::sym("x", states->nx) somewhere
                casadi::SX u; // this needs to be initialized to casadi::SX::sym("u", states->nu) somewhere
                casadi::SX t; // this needs to be initialized to casadi::SX::sym("t") somewhere
                int num_knots;
                int collocation_points_per_knot; // TODO get rid of this

                double max_footstep_offset_height;
                double corrector_kp;
                double following_leeway = 0;
            };

            template <class ProblemData>
            class VelocityConstraintBuilder : Constraint<ProblemData>
            {

            public:
                VelocityConstraintBuilder() : ConstraintBuilder() {}

                /**
                 *
                 * @brief Generate flags for each knot point. We set it to all ones, applicable at each knot.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "VelocityConstraintProblemData" NAMED "velocity_constraint_problem_data"
                 * @param apply_at
                 */
                void CreateApplyAt(const ProblemData &problem_data, int knot_index, Eigen::VectorXi &apply_at) const override
                {
                    uint num_points = problem_data.friction_cone_problem_data.collocation_points_per_knot;
                    apply_at = Eigen::VectorXi::Constant(num_points, 1);
                }

                void BuildConstraint(const ProblemData &problem_data, int knot_index, ConstraintData &constraint_data) const override;
            };

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::BuildConstraint(const ProblemData &problem_data, int knot_index, ConstraintData &constraint_data) const
            {
                std::vector<casadi::Function> G_vec;
                std::vector<casadi::Function> upper_bound_vec;
                std::vector<casadi::Function> lower_bound_vec;

                CreateApplyAt(problem_data, knot_index, constraint_data.apply_at);

                casadi::SX t("t", 1);

                ContactSequence::CONTACT_SEQUENCE_ERROR error;

                int phase_index_at_knot = problem_data.contact_sequence->getPhaseIndexAtKnot(knot_index, error);
                auto mode = problem_data.contact_sequence->getPhase(phase_index_at_knot).mode;

                for (int ee_index = 0; ee_index < problem_data.robot_end_effectors.size(); ee_index++)
                {
                    auto ee = (problem_data.robot_end_effectors[ee_index].second);
                    if (!mode[*ee])
                    {
                        double liftoff_time = 0;
                        double touchdown_time = problem_data.contact_sequence->dt();
                        // When did the EE break contact?
                        for (int i = phase_index_at_knot; i >= 0; i--)
                        {
                            // If the ee is in contact, we have found the phase where liftoff occurred.
                            if (problem_data.contact_sequence->getPhase(i).mode[*ee])
                            {
                                int liftoff_index = i + 1;
                                problem_data.contact_sequence->getTimeAtPhase(liftoff_index, liftoff_time, error);
                                break;
                            }
                        }

                        for (int i = phase_index_at_knot; i < problem_data.contact_sequence->numPhases() - 1; i++)
                        {
                            // If the ee is in contact, we have found the phase where touchdown occurred.
                            if (problem_data.contact_sequence->getPhase(i).mode[*ee])
                            {
                                int touchdown_index = i;
                                problem_data.contact_sequence->getTimeAtPhase(touchdown_index, touchdown_time, error);
                                break;
                                // if touchdown is not found, we will assume that the ee is in contact at the end of the horizon.
                                //  This is an odd and limiting assumption. Better behavior should be created.
                            }
                        }
                        // Add a baumgarte corrector height_velocity = (kp * (position(x) - desired_position) + desired_velocity)

                        // velocity_z - kp * (position_z(state)) = desired_velocity_z - kp * desired_position_z

                        // velocity_z - kp * (position_z(state)) < desired_velocity_z - kp * desired_position_z + leeway
                        // velocity_z - kp * (position_z(state)) > desired_velocity_z - kp * desired_position_z - leeway

                        casadi::Function t_in_range = casadi::Function("t_in_range", {t}, {(t - liftoff_time) / (touchdown_time - liftoff_time)});

                        auto footstep_height_function = createFootstepHeightFunction(t_in_range, problem_data.footstep_trajectory_generator->getFootstepDefinition(*ee));
                        casadi::Function desired_height = footstep_height_function[0];
                        casadi::Function desired_velocity = footstep_height_function[1];

                        casadi::Function vel_bound = desired_velocity - problem_data.corrector_kp * desired_height;

                        auto &frame_omf_data = problem_data.velocity_constraint_problem_data->ad_data->oMf[ee->frame_id];

                        casadi::Function foot_pos = frame_omf_data.translation();
                        casadi::Function foot_vel = frame_omf_data.velocity().linear();

                        casadi::Function vel_constraint_G = foot_vel - problem_data.corrector_kp * foot_pos;

                        G_vec.push_back(vel_constraint_G);
                        lower_bound_vec.push_back(vel_bound - problem_data.following_leeway);
                        upper_bound_vec.push_back(vel_bound + problem_data.following_leeway);
                    }
                    else
                    {
                        // add velocity = 0 as a constraint

                        casadi::SX foot_vel = frame_omf_data.velocity().linear();

                        G_vec.push_back(vel_constraint_G);
                        lower_bound_vec.push_back(casadi::SX::zeros(3, 1));
                        upper_bound_vec.push_back(casadi::SX::zeros(3, 1));
                    }
                }

                constraint_data.G = vertcat(G_vec);
                constraint_data.upper_bound = vertcat(upper_bound_vec);
                constraint_data.lower_bound = vertcat(lower_bound_vec);
            }

        }
    }
}