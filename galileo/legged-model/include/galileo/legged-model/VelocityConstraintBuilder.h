#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"
#include "galileo/legged-model/LeggedRobotStates.h"

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

                std::shared_ptr<FootstepTrajectoryGenerator> footstep_trajectory_generator;
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

                /**
                 * @brief Generate bounds for a vector of points.
                 *
                 * For both approximations, each value is less than 0, with no lower bound.
                 *t each knot.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "VelocityConstraintProblemData" NAMED "velocity_constraint_problem_data"
                 * @param upper_bound
                 * @param lower_bound
                 */
                void CreateBounds(const ProblemData &problem_data, int knot_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const override;

                /**
                 * @brief Generate a function to evaluate each point
                 *
                 * t each knot.
                 *
                 * @param problem_data MUST CONTAIN AN INSTANCE OF "VelocityConstraintProblemData" NAMED "velocity_constraint_problem_data"
                 * @param G
                 */
                void CreateFunction(const ProblemData &problem_data, int knot_index, casadi::Function &G) const;

            private:
                /**
                 * @brief Get the mode at a given knot index.
                 *
                 * @param problem_data
                 * @param knot_index
                 * @return ContactMode
                 */
                const contact::ContactMode &getModeAtKnot(const ProblemData &problem_data, int knot_index) const;

                /**
                 * @brief Get the desired EE velocity for all collocation points in a knot.
                 */
                std::vector<Eigen::Vector2d> getDesiredEEHeightState(const ProblemData &problem_data, std::shared_ptr<contact::EndEffector> ee, int contact_broken_index, int contact_made_index);

                int getNumberOfConstraintsPerEE(const ProblemData &problem_data) const
                {
                    // Assuming ONE constraint on velocity per collocation point. If there are constraints on all three dimensions, then this needs to be changed.
                    return problem_data.collocation_points_per_knot;
                }
            };

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::CreateFunction(const ProblemData &problem_data, int knot_index, casadi::Function &G) const;
            {
                auto mode = getModeAtKnot(problem_data, knot_index);

                for (int ee_index = 0; ee_index < problem_data.robot_end_effectors.size(); ee_index++)
                {
                    auto ee = (problem_data.robot_end_effectors[ee_index].second);
                    if (!mode[*ee])
                    {
                        // When did the EE break contact?
                        int contact_broken_index = 0;
                        int contact_made_index = problem_data.num_knots - 1;
                        for (int i = knot_index; i > 0; i--)
                        {
                            if ((getModeAtKnot(problem_data, i)[*ee]))
                            {
                                contact_broken_index = i + 1;
                                break;
                            }
                        }

                        // When will the EE make contact again?
                        for (int i = knot_index; i < problem_data.num_knots; i++)
                        {
                            if ((getModeAtKnot(problem_data, i)[*ee]))
                            {
                                contact_made_index = i - 1;
                                break;
                            }
                        }
                    }

                    std::vector<Eigen::Vector2d> desired_ee_height_state = getDesiredEEHeightState(problem_data, ee, contact_broken_index, contact_made_index);

                    // Add a baumgarte corrector height_velocity = (kp * (position(x) - desired_position) + desired_velocity)
                }
            }

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::CreateBounds(const ProblemData &problem_data, int knot_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
            {
                int num_constraints_per_ee = getNumberOfConstraintsPerEE(problem_data);
                int num_constraints = num_constraints_per_ee * problem_data.robot_end_effectors.size();

                casadi::SX upper_bound_sx(num_constraints);
                casadi::SX lower_bound_sx(num_constraints);

                upper_bound_sx = casadi::SX::zeros(num_constraints);
                lower_bound_sx = casadi::SX::zeros(num_constraints);

                for (int ee_index = 0; ee_index < problem_data.robot_end_effectors.size(); ee_index++)
                {
                    auto ee = (problem_data.robot_end_effectors[ee_index].second);
                    if (!mode[*ee])
                    {
                        // We add some leeway. The velocity can be slightly above or below the desired velocity. If this is used, then there should be a constraint on the initial and final height _position_.

                        // End Effector i has constraints (ee_index * num_constraints_per_ee)
                        int ee_constraint_index_start = ee_index * num_constraints_per_ee;
                        int ee_constraint_index_end = ee_constraint_index_start + num_constraints_per_ee;

                        // upper_bound_sx.block(ee_constraint_index_start, ee_constraint_index_end) = zeros(num_constraints_per_ee) + problem_data.following_leeway;
                        // upper_bound_sx.block(ee_constraint_index_start, ee_constraint_index_end) = zeros(num_constraints_per_ee) - problem_data.following_leeway;
                    }
                }

                upper_bound = casadi::Function("upper_bound", {problem_data.x, problem_data.u, problem_data.t}, {upper_bound_sx});
                lower_bound = casadi::Function("lower_bound", {problem_data.x, problem_data.u, problem_data.t}, {lower_bound_sx});
            }

            template <class ProblemData>
            std::vector<Eigen::Vector2d> VelocityConstraintBuilder<ProblemData>::getDesiredEEHeightState(const ProblemData &problem_data, std::shared_ptr<contact::EndEffector> ee, int contact_broken_index, int contact_made_index)
            {
                // How many collocation points are there between the contact breaking and making?
                int broken_contact_span =
                    (contact_made_index - contact_broken_index) * problem_data.collocation_points_per_knot;

                // How many collocation points are there between the current knot and the contact breaking?
                int broken_to_current_span =
                    (knot_index - contact_broken_index) * problem_data.collocation_points_per_knot;

                FootstepTrajectoryGenerator::FootstepDefinition fs_def;

                // The liftoff and touchdown surfaces for the end effector
                environment::SurfaceID pre_surface_id = getModeAtKnot(problem_data, contact_broken_index).getSurfaceID(ee);
                environment::SurfaceID post_surface_id = getModeAtKnot(problem_data, contact_made_index).getSurfaceID(ee);

                // The liftoff and touchdown heights for the end effector
                if (pre_surface_id == environment::NO_SURFACE)
                {
                    // Need to define behavior for this case.
                    fs_def.h_start = 0;
                }
                else
                {
                    fs_def.h_start = (*problem_data.environment_surfaces)[pre_surface_id].height;
                }

                if (post_surface_id == environment::NO_SURFACE)
                {
                    // Need to define behavior for this case.
                    fs_def.h_end = 0;
                }
                else
                {
                    fs_def.h_end = (*problem_data.environment_surfaces)[post_surface_id].height;
                }

                fs_def.h_max = problem_data.max_footstep_offset_height + std::min(fs_def.h_start, fs_def.h_end);

                asset(fs_def.h_max > std::max(fs_def.h_start, fs_def.h_end));

                std::vector<Eigen::Vector2d> desired_ee_height_state(broken_contact_span);
                for (int i = 0; i < broken_contact_span; i++)
                {
                    desired_ee_height_state[i] = problem_data.footstep_trajectory_generator->calcFootstepZ(broken_to_current_span + i, broken_contact_span, fs_def);
                }

                return desired_ee_height_state;
            }

        }
    }
}