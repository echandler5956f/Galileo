#pragma once

#include "galileo/opt/Constraint.h"
#include "galileo/legged-model/ContactSequence.h"
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {

            struct FootstepDefinition
            {
                double h_offset_P1;
                double h_offset_P2;

                // R s.t. R * z_hat = normal of surface
                Eigen::MatrixXd R_P1;
                Eigen::MatrixXd R_P2;

                // A value between 0 and 1 that determines the window of time on which the footstep has velocitty along basis n.
                // We consider the total spatial velocity to be some [n1 n2] * [v1; v2].
                // At t = h_1_window_duration, the velocity v1 is 0.  at t = 1 - h2_window_duration, the velocity v2 is 0.
                double h1_window_duration = 0.5;
                double h2_window_duration = 0.5;

                // The number of standard deviations we consider in the window of the bell curve.
                // Effectively, this defines the steepness fo the bell curve.
                double window_sigma = 3.3; // 3.3

                // How tightly thhe bell curve trajectory is followed. The higher the value, the more tightly the trajectory is followed.
                double sigmoid_scaling = 15.0;

                double liftoff_time;
                double touchdown_time;
            };

            casadi::Function getMappedBellCurve(casadi::SX &t, double window_sigma)
            {
                // The bell curve is defined as a function of t and mu, where t is the time from 0 to 1.
                casadi::Function bell_curve = casadi::Function("bell_curve", casadi::SXVector{t}, casadi::SXVector{exp(-pow((t), 2) / (2))});
                // The bell curve mapped such that the mean is at 0.5, mu - window_sigma is at 0 and mu + window_sigma is at 1. Standard deviation is assumed to be 1
                casadi::Function mapped_bell_curve = casadi::Function("mapped_bell_curve", {t},
                                                                      {2 * window_sigma * bell_curve(casadi::SXVector{(t - 0.5) * (2 * window_sigma)}).at(0) / (sqrt(2 * M_PI))});
                return mapped_bell_curve;
            }

            casadi::Function getSigmoid(casadi::SX &t, double sigmoid_scaling)
            {
                casadi::Function sigmoid = casadi::Function("sigmoid", casadi::SXVector{t}, casadi::SXVector{1 / (1 + exp(-sigmoid_scaling * (t - 0.5)))});
                casadi::Function fixed_sigmoid = casadi::Function("fixed_sigmoid", casadi::SXVector{t},
                                                                  casadi::SXVector{
                                                                      sigmoid(casadi::SXVector{t}).at(0) - sigmoid(casadi::SXVector{0}).at(0)});
                return sigmoid;
            }

            /**
             * @brief A struct for holding the velocity constraint problem data.
             */
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

                // The ideal height a footstep should be off a surface
                double ideal_offset_height;

                double max_following_leeway_planar = 0.5;
                double min_following_leeway_planar = 1e-5;

                double max_following_leeway_normal = 0.5;
                double min_following_leeway_normal = 1e-5;

                // How tightly thhe bell curve trajectory is followed. The higher the value, the more tightly the trajectory is followed.
                double sigmoid_scaling = 15.0;
            };

            /**
             * @brief A class for building the velocity constraint.
             */
            template <class ProblemData>
            class VelocityConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
            {

            public:
                VelocityConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

                void buildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data);

                void createBounds(const ProblemData &problem_data, int phase_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
                {
                }

                void createFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
                {
                }

                /**
                 * @brief The velocity of the footstep in the direction of the surface normals.
                 *
                 * @param t The time of the trajectory (casadi symbolic variable)
                 * @param FS_def The footstep definition.
                 *
                 * @return The desired velocity of the footstep in the direction of the surface normals.
                 *
                 */
                void footstepVelocityFunction(casadi::SX &t, const FootstepDefinition &FS_def, casadi::Function &h1_dot_desired, casadi::Function &h2_dot_desired)
                {
                    // Volocity in world frame = h_dot_1 * n1 + h_dot_2 * n2 where h_dot_n gives the velocity of the footstep in the direction of n

                    // We will define the velocity as a bell curve.
                    // The bell curve will be defined as a function of t, where t is the time from 0 to 1.
                    // The mean of the bell curve represents the time when velocity is highest.
                    // The standard deviation defines the steepness of the bell curve.

                    double window_sigma = FS_def.window_sigma;

                    casadi::Function mapped_bell_curve = getMappedBellCurve(t, window_sigma);

                    h1_dot_desired = casadi::Function("h1_dot_desired", {t}, {FS_def.h_offset_P1 * mapped_bell_curve(casadi::SXVector{t / FS_def.h1_window_duration}).at(0)});
                    h2_dot_desired = casadi::Function("h2_dot_desired", {t}, {FS_def.h_offset_P2 * mapped_bell_curve(casadi::SXVector{(1 - t) / FS_def.h2_window_duration}).at(0)});
                }

                casadi::SX getFootstepVelocityInWorldFrame(const ProblemData &problem_data, pinocchio::FrameIndex frame_id) const
                {
                    Eigen::Matrix<galileo::opt::ADScalar, 6, 1, 0> foot_vel = pinocchio::getFrameVelocity(*(problem_data.velocity_constraint_problem_data.ad_model),
                                                                                                          *(problem_data.velocity_constraint_problem_data.ad_data),
                                                                                                          frame_id,
                                                                                                          pinocchio::LOCAL)
                                                                                  .toVector();

                    casadi::SX cfoot_vel = casadi::SX(casadi::Sparsity::dense(foot_vel.rows(), 1));
                    pinocchio::casadi::copy(foot_vel, cfoot_vel);

                    return cfoot_vel;
                }

                casadi::SX getFootstepVelocityInSurfaceFrame(casadi::SX &v, const Eigen::MatrixXd &R_n) const
                {
                    casadi::SX R_n_inv_casadi = casadi::SX::zeros(3, 3);
                    pinocchio::casadi::copy(R_n.transpose(), R_n_inv_casadi);

                    casadi::SX v_in_surface_frame = casadi::SX::mtimes(R_n_inv_casadi, v);
                    return v_in_surface_frame;
                }

                FootstepDefinition buildFootstepDefinition(const ProblemData &problem_data, int phase_index, std::shared_ptr<galileo::legged::contact::EndEffector> ee_ptr) const;
            };

            template <class ProblemData>
            FootstepDefinition VelocityConstraintBuilder<ProblemData>::buildFootstepDefinition(const ProblemData &problem_data, int phase_index, std::shared_ptr<galileo::legged::contact::EndEffector> ee_ptr) const
            {

                double liftoff_time = 0;
                double touchdown_time = problem_data.velocity_constraint_problem_data.contact_sequence->getDT();
                FootstepDefinition footstep_definition;

                galileo::legged::environment::SurfaceID liftoff_surface_ID = 0;
                galileo::legged::environment::SurfaceID touchdown_surface_ID = 0;

                contact::ContactSequence::CONTACT_SEQUENCE_ERROR error;

                // When did the EE break contact?
                for (int i = phase_index; i >= 0; i--)
                {
                    contact::ContactMode mode_i = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode;
                    // If the ee is in contact, we have found the phase where liftoff occurred.
                    if (mode_i.at(*ee_ptr))
                    {
                        int liftoff_index = i + 1;
                        problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(liftoff_index, liftoff_time, error);
                        liftoff_surface_ID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee_ptr);
                        break;
                        // if liftoff is not found, we will assume that the ee is in contact with SURFACE 0 at the start of the horizon.
                        //  This is an odd and limiting assumption. Better behavior should be created.
                    }
                }

                // When does it make contact again?
                for (int i = phase_index; i < problem_data.velocity_constraint_problem_data.contact_sequence->getNumPhases(); i++)
                {
                    contact::ContactMode mode_i = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode;
                    // If the ee is in contact, we have found the phase where touchdown occurred.
                    if (mode_i.at(*ee_ptr))
                    {
                        int touchdown_index = i;
                        problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(touchdown_index, touchdown_time, error);
                        touchdown_surface_ID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee_ptr);
                        break;
                        // if liftoff is not found, we will assume that the ee is in contact with SURFACE 0 at the start of the horizon.
                        //  This is an odd and limiting assumption. Better behavior should be created.
                    }
                }
                footstep_definition.h_offset_P1 = problem_data.velocity_constraint_problem_data.ideal_offset_height;
                footstep_definition.h_offset_P2 = problem_data.velocity_constraint_problem_data.ideal_offset_height;

                footstep_definition.R_P1 = problem_data.velocity_constraint_problem_data.environment_surfaces->getSurfaceFromID(liftoff_surface_ID).Rotation();
                footstep_definition.R_P2 = problem_data.velocity_constraint_problem_data.environment_surfaces->getSurfaceFromID(touchdown_surface_ID).Rotation();

                footstep_definition.liftoff_time = liftoff_time;
                footstep_definition.touchdown_time = touchdown_time;

                footstep_definition.sigmoid_scaling = problem_data.velocity_constraint_problem_data.sigmoid_scaling;

                return footstep_definition;
            }

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::buildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data)
            {
                casadi::SXVector G_vec;
                casadi::SXVector upper_bound_vec;
                casadi::SXVector lower_bound_vec;

                casadi::SX t = problem_data.velocity_constraint_problem_data.t;

                contact::ContactMode mode = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(phase_index).mode;

                for (auto ee : problem_data.velocity_constraint_problem_data.robot_end_effectors)
                {
                    if (!mode[(*ee.second)])
                    {
                        FootstepDefinition footstep_definition = buildFootstepDefinition(problem_data, phase_index, ee.second);

                        /**
                         *
                         *
                         * the velocity _in the frame of the liftoff and touchdown surfaces_ is constrained instead.
                         * The normal velocity from the liftoff surface is constrained such that it approximately follows a bell curve.
                         * As time goes on, however, the constraint is relaxed, and the velocity normal is allowed to be greater.
                         * Similarly, the normal velocity from the touchdown surface is constrained such that it approximately follows a bell curve.
                         * At t = 0, however, the constraint is less restrictive.
                         *
                         * In essence, the velocity at t = 0 approximately is normal to the liftoff surface and at t = 1, approximately normal to the touchdown surface.
                         * In the middle, it loosely follows an "ideal" trajectory.
                         */
                        casadi::SX cfoot_vel = getFootstepVelocityInWorldFrame(problem_data, ee.first)(casadi::Slice(0, 3), 0);
                        casadi::SX cfoot_vel_in_liftoff_surface = getFootstepVelocityInSurfaceFrame(cfoot_vel, footstep_definition.R_P1);
                        casadi::SX cfoot_vel_in_touchdown_surface = getFootstepVelocityInSurfaceFrame(cfoot_vel, -footstep_definition.R_P2);

                        // The velocity of the footstep in the direction of each of the surface normals
                        // At time 0, the velocity should be approximately normal to the liftoff surface
                        // At time 1, the velocity should be approximately normal to the touchdown surface
                        casadi::Function G_foot_vel = casadi::Function("G_foot_vel", casadi::SXVector{problem_data.velocity_constraint_problem_data.x, problem_data.velocity_constraint_problem_data.u},
                                                                       casadi::SXVector{vertcat(casadi::SXVector{cfoot_vel_in_liftoff_surface, cfoot_vel_in_touchdown_surface})});

                        G_vec.push_back(G_foot_vel(casadi::SXVector{problem_data.velocity_constraint_problem_data.x, problem_data.velocity_constraint_problem_data.u}).at(0));

                        // Building the bounds. We want to approximately follow (\dot h1_desired) for a fraction of the trajectory, then (\dot h2_desired).

                        // Map t to between 0 and 1. 0 is the time of liftoff, 1 is the time of touchdown.
                        casadi::SX mapped_t = (t - footstep_definition.liftoff_time) / (footstep_definition.touchdown_time - footstep_definition.liftoff_time);

                        // Desired velocity normal to the liftoff surface as a function of time
                        casadi::Function desired_h1_dot;
                        casadi::Function desired_h2_dot;

                        footstepVelocityFunction(t, footstep_definition, desired_h1_dot, desired_h2_dot);

                        casadi::SX ell_max_planar = problem_data.velocity_constraint_problem_data.max_following_leeway_planar;
                        casadi::SX ell_min_planar = problem_data.velocity_constraint_problem_data.min_following_leeway_planar;
                        casadi::SX ell_slope_planar = (ell_max_planar - ell_min_planar);

                        casadi::SX ell_max_normal = problem_data.velocity_constraint_problem_data.max_following_leeway_normal;
                        casadi::SX ell_min_normal = problem_data.velocity_constraint_problem_data.min_following_leeway_normal;
                        casadi::SX ell_slope_normal = (ell_max_normal - ell_min_normal);

                        // Linear interpolation between admissible error of velocity at t = [0, 1]
                        casadi::Function linear_error_interpolation_planar = casadi::Function("linear_error_interpolation_planar", {t}, {casadi::SX(ell_slope_planar * t + ell_min_planar)});
                        casadi::Function linear_error_interpolation_normal = casadi::Function("linear_error_interpolation_normal", {t}, {casadi::SX(ell_slope_normal * t + ell_min_normal)});

                        // The admissible error of velocity parallel to the liftoff surface. As t approaches 0, this should approach ell_min
                        casadi::SX admissible_error_h1_parallel = linear_error_interpolation_planar(t).at(0);
                        // The admissible error of velocity parallel to the touchdown surface. As t approaches 1, this should approach ell_min
                        casadi::SX admissible_error_h2_parallel = linear_error_interpolation_planar(1 - t).at(0);

                        // The admissible error of velocity normal to the liftoff surface. As t approaches 0, this should approach ell_min
                        casadi::SX upper_admissible_error_h1_normal = linear_error_interpolation_normal(t).at(0);
                        // The admissible error of velocity normal to the touchdown surface. As t approaches 1, this should approach ell_min
                        casadi::SX upper_admissible_error_h2_normal = linear_error_interpolation_normal(1 - t).at(0);

                        // The lower bound of the velocity normal to the surfaces. We want this to evolve slower, so that the velocity is "encouraged" to be positive
                        //  This is the "magnitude" of the offset from h_dot_desired. The actual bound is h_dot_desired - lower_admissible_error_h_normal
                        casadi::Function quadratic_error_interpolation = casadi::Function("quadratic_error_interpolation", {t}, {casadi::SX(ell_slope_normal * pow(2 * t, 2) + ell_min_normal)});

                        casadi::Function sigmoid = getSigmoid(t, footstep_definition.sigmoid_scaling);
                        casadi::SX max_error_offset = desired_h1_dot(casadi::SXVector{footstep_definition.h1_window_duration / 2}).at(0) + 0.2;
                        casadi::SX lower_admissible_error_h1_normal = max_error_offset * sigmoid(t).at(0) + ell_min_normal;
                        casadi::SX lower_admissible_error_h2_normal = max_error_offset * sigmoid(1 - t).at(0) + ell_min_normal;
                        //  casadi::SX lower_admissible_error_h1_normal = upper_admissible_error_h1_normal;
                        // casadi::SX lower_admissible_error_h2_normal = upper_admissible_error_h2_normal;

                        casadi::Function upper_bound = casadi::Function("upper_bound", casadi::SXVector{t},
                                                                        casadi::SXVector{vertcat(casadi::SXVector{
                                                                            admissible_error_h1_parallel,
                                                                            admissible_error_h1_parallel,
                                                                            desired_h1_dot(t).at(0) + upper_admissible_error_h1_normal,
                                                                            admissible_error_h2_parallel,
                                                                            admissible_error_h2_parallel,
                                                                            desired_h2_dot(t).at(0) + upper_admissible_error_h2_normal})});

                        casadi::Function lower_bound = casadi::Function("lower_bound", casadi::SXVector{t},
                                                                        casadi::SXVector{vertcat(casadi::SXVector{
                                                                            -admissible_error_h1_parallel,
                                                                            -admissible_error_h1_parallel,
                                                                            desired_h1_dot(t).at(0) - lower_admissible_error_h1_normal,
                                                                            -admissible_error_h2_parallel,
                                                                            -admissible_error_h2_parallel,
                                                                            desired_h2_dot(t).at(0) - lower_admissible_error_h2_normal})});

                        lower_bound_vec.push_back(lower_bound(mapped_t).at(0));
                        upper_bound_vec.push_back(upper_bound(mapped_t).at(0));
                    }
                    else
                    {
                        casadi::SX cfoot_vel = getFootstepVelocityInWorldFrame(problem_data, ee.first);
                        if (ee.second->is_6d)
                        {
                            G_vec.push_back(cfoot_vel);
                            lower_bound_vec.push_back(casadi::SX::zeros(cfoot_vel.rows(), 1));
                            upper_bound_vec.push_back(casadi::SX::zeros(cfoot_vel.rows(), 1));
                        }
                        else
                        {
                            G_vec.push_back(cfoot_vel(casadi::Slice(0, 3), 0));
                            lower_bound_vec.push_back(casadi::SX::zeros(3, 1));
                            upper_bound_vec.push_back(casadi::SX::zeros(3, 1));
                        }
                    }
                }
                constraint_data.G = casadi::Function("G_Velocity", casadi::SXVector{problem_data.velocity_constraint_problem_data.x, problem_data.velocity_constraint_problem_data.u}, casadi::SXVector{vertcat(G_vec)});
                constraint_data.upper_bound = casadi::Function("upper_bound", casadi::SXVector{t}, casadi::SXVector{vertcat(upper_bound_vec)});
                constraint_data.lower_bound = casadi::Function("lower_bound", casadi::SXVector{t}, casadi::SXVector{vertcat(lower_bound_vec)});
                std::cout << "G_Velocity: " << constraint_data.G << std::endl;
                std::cout << "upper_bound: " << constraint_data.upper_bound << std::endl;
                std::cout << "lower_bound: " << constraint_data.lower_bound << std::endl;
            }

        }
    }
}