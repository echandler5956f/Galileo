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
#include <optional>

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {

            /**
             * @brief Struct for defining the footstep.
             *
             */
            struct FootstepDefinition
            {
                /**
                 * @brief Liftoff time of the footstep.
                 *
                 */
                double liftoff_time;

                /**
                 * @brief Touchdown time of the footstep.
                 *
                 */
                double touchdown_time;

                /**
                 * @brief Velocity at the start of the footstep.
                 *
                 */
                double footstep_vel_start;

                /**
                 * @brief Velocity at the end of the footstep.
                 *
                 */
                double footstep_vel_end;

                /**
                 * @brief Start height of the footstep.
                 *
                 */
                double h_start;

                /**
                 * @brief End height of the footstep.
                 *
                 */
                double h_end;

                /**
                 * @brief Maximum height of the footstep.
                 *
                 */
                double h_max;
            };

            casadi::SX createFootstepHeightFunction(casadi::SX t, const FootstepDefinition FS_def, double scaled)
            {
                double takeoff_position = FS_def.h_start;
                double takeoff_velocity = FS_def.footstep_vel_start * scaled;
                double midheight = FS_def.h_max * scaled / (0.5 * (FS_def.touchdown_time - FS_def.liftoff_time));
                double mid_velocity = 0. * scaled;
                double touchdown_position = FS_def.h_end;
                double touchdown_velocity = FS_def.footstep_vel_end;

                casadi::SX terms = vertcat(casadi::SXVector{1, t, pow(t, 2), pow(t, 3)});
                casadi::SX A = vertcat(casadi::SXVector{horzcat(casadi::SXVector{1, 0, 0, 0}),
                                                        horzcat(casadi::SXVector{0, 1, 0, 0}),
                                                        horzcat(casadi::SXVector{1, 1, 1, 1}),
                                                        horzcat(casadi::SXVector{0, 1, 2, 3})});
                casadi::SX b_takeoff = vertcat(casadi::SXVector{takeoff_position, takeoff_velocity,
                                                                midheight, mid_velocity});
                casadi::SX x_takeoff = casadi::SX::mtimes(casadi::SX::inv(A), b_takeoff);
                casadi::SX p_takeoff = casadi::SX::mtimes(terms.T(), x_takeoff);
                casadi::Function p_takeoff_func = casadi::Function("p1_func", {t}, {p_takeoff});
                casadi::SX dp_takeoff = casadi::SX::jacobian(p_takeoff, t);
                casadi::Function dp_takeoff_func = casadi::Function("dp1_func", {t}, {dp_takeoff});

                casadi::SX b_touchdown = vertcat(casadi::SXVector{midheight, mid_velocity,
                                                                  touchdown_position, touchdown_velocity});
                casadi::SX x_touchdown = casadi::SX::mtimes(casadi::SX::inv(A), b_touchdown);
                casadi::SX p_touchdown = casadi::SX::mtimes(terms.T(), x_touchdown);
                casadi::Function p_touchdown_func = casadi::Function("p2_func", {t}, {p_touchdown});
                casadi::SX dp_touchdown = casadi::SX::jacobian(p_touchdown, t);
                casadi::Function dp_touchdown_func = casadi::Function("dp2_func", {t}, {dp_touchdown});

                casadi::SX piecewise_poly = casadi::SX::if_else(t < 0.5, dp_takeoff_func(casadi::SXVector{t / 0.5}).at(0), dp_touchdown_func(casadi::SXVector{(t - 0.5) / (1 - 0.5)}).at(0));
                return piecewise_poly;
            }

            /**
             * @brief A struct for holding the velocity constraint problem data.
             */
            struct VelocityConstraintProblemData
            {
                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<legged::LeggedRobotStates> states;
                std::shared_ptr<legged::ADModel> ad_model;
                std::shared_ptr<legged::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;

                // The height a footstep should be off a surface
                double ideal_offset_height = 0.08;
                double footstep_height_scaling = 0.15;
                double max_following_leeway_planar = 50;
                double min_following_leeway_planar = 1e-8;
                double footstep_vel_start = 0;
                double footstep_vel_end = 0;

                std::optional<std::map<std::string, double>> footstep_liftoffs_before_horizon = std::nullopt;
                std::optional<std::map<std::string, double>> footstep_touchdowns_after_horizon = std::nullopt;

                std::optional<double> ideal_footstep_duration = std::nullopt;
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

                casadi::SX getFootstepVelocity(const ProblemData &problem_data, pinocchio::FrameIndex frame_id) const
                {
                    Eigen::Matrix<legged::ADScalar, 6, 1, 0> foot_vel = pinocchio::getFrameVelocity(*(problem_data.velocity_constraint_problem_data.ad_model),
                                                                                                    *(problem_data.velocity_constraint_problem_data.ad_data),
                                                                                                    frame_id,
                                                                                                    pinocchio::LOCAL_WORLD_ALIGNED)
                                                                            .toVector();

                    casadi::SX cfoot_vel = casadi::SX(casadi::Sparsity::dense(foot_vel.rows(), foot_vel.cols()));
                    pinocchio::casadi::copy(foot_vel, cfoot_vel);

                    return cfoot_vel;
                }
            };

            template <class ProblemData>
            void VelocityConstraintBuilder<ProblemData>::buildConstraint(const ProblemData &problem_data, int phase_index, opt::ConstraintData &constraint_data)
            {
                casadi::SXVector G_vec;
                casadi::SXVector lower_bound_vec;
                casadi::SXVector upper_bound_vec;

                casadi::SX t = problem_data.velocity_constraint_problem_data.t;
                contact::ContactMode mode = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(phase_index).mode;

                auto liftoffs_before_horizon = problem_data.velocity_constraint_problem_data.footstep_liftoffs_before_horizon;
                auto touchdowns_after_horizon = problem_data.velocity_constraint_problem_data.footstep_touchdowns_after_horizon;
                for (auto ee : problem_data.velocity_constraint_problem_data.robot_end_effectors)
                {
                    if (!mode[(*ee.second)])
                    {
                        double liftoff_time = 0;
                        double touchdown_time = problem_data.velocity_constraint_problem_data.contact_sequence->getDT();
                        FootstepDefinition footstep_definition;
                        footstep_definition.h_start = 0;
                        footstep_definition.h_end = 0;
                        footstep_definition.footstep_vel_start = problem_data.velocity_constraint_problem_data.footstep_vel_start;
                        footstep_definition.footstep_vel_end = problem_data.velocity_constraint_problem_data.footstep_vel_end;

                        galileo::legged::environment::SurfaceID liftoff_surface_ID = 0;
                        galileo::legged::environment::SurfaceID touchdown_surface_ID = 0;

                        contact::ContactSequence::CONTACT_SEQUENCE_ERROR error;

                        // When did the EE break contact?
                        for (int i = phase_index; i >= 0; i--)
                        {
                            contact::ContactMode mode_i = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode;
                            // If the ee is in contact, we have found the phase where liftoff occurred.
                            if (mode_i.at(*ee.second))
                            {
                                int liftoff_index = i + 1;
                                problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(liftoff_index, liftoff_time, error);
                                liftoff_surface_ID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee.second);
                                break;
                            }
                        }


                        // When does it make contact again?
                        for (int i = phase_index; i < problem_data.velocity_constraint_problem_data.contact_sequence->getNumPhases(); i++)
                        {
                            contact::ContactMode mode_i = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode;
                            // If the ee is in contact, we have found the phase where touchdown occurred.
                            if (mode_i.at(*ee.second))
                            {
                                int touchdown_index = i;
                                problem_data.velocity_constraint_problem_data.contact_sequence->getTimeAtPhase(touchdown_index, touchdown_time, error);
                                touchdown_surface_ID = problem_data.velocity_constraint_problem_data.contact_sequence->getPhase(i).mode.getSurfaceID(*ee.second);
                                break;
                            }
                        }

                        if(liftoff_time == 0){
                            // The EE is not in contact at the start of the horizon
                            if(liftoffs_before_horizon.has_value() && liftoffs_before_horizon.value().count(ee.second->frame_name) > 0){
                                liftoff_time =  liftoffs_before_horizon.value()[ee.second->frame_name];
                            } else {
                                if(touchdown_time < problem_data.velocity_constraint_problem_data.contact_sequence->getDT()){
                                    if(problem_data.velocity_constraint_problem_data.ideal_footstep_duration.has_value()){
                                        liftoff_time = touchdown_time - problem_data.velocity_constraint_problem_data.ideal_footstep_duration.value();
                                    } else {
                                        throw std::runtime_error("The ideal footstep liftoff time for " + ee.second->frame_name + " is not defined!");
                                    }
                                } else{
                                    // The end effector is in flight for the entire horizon! The velocity constraint cannot be applied
                                    
                                }
                            }
                        
                            if(liftoff_time > 0){
                                // The liftoff is within the horizon, but the mode at t=0 is not in contact
                                // this is an error
                                throw std::runtime_error("The ideal footstep liftoff time for " + ee.second->frame_name + " is defined, but it is within the planning horizon and the EE is not in contact at that time!");
                            }
                        }

                        if(touchdown_time == problem_data.velocity_constraint_problem_data.contact_sequence->getDT()){
                            // The EE is not in contact at the end of the horizon
                            if(touchdowns_after_horizon.has_value() && touchdowns_after_horizon.value().count(ee.second->frame_name) > 0){
                                touchdown_time =  touchdowns_after_horizon.value()[ee.second->frame_name];
                            } else {
                                if(liftoff_time > 0){
                                    if(problem_data.velocity_constraint_problem_data.ideal_footstep_duration.has_value()){
                                        touchdown_time = liftoff_time + problem_data.velocity_constraint_problem_data.ideal_footstep_duration.value();
                                    } else {
                                        throw std::runtime_error("The ideal footstep touchdown time for " + ee.second->frame_name + " is not defined!");
                                    }
                                } else{
                                    // The end effector is in flight for the entire horizon! The velocity constraint cannot be applied
                                    
                                }
                            }
                        
                            if(touchdown_time < problem_data.velocity_constraint_problem_data.contact_sequence->getDT()){
                                // The touchdown is within the horizon, but the mode at t=DT is not in contact
                                // this is an error
                                throw std::runtime_error("The ideal footstep touchdown time for " + ee.second->frame_name + " is defined, but it is within the planning horizon and the EE is not in contact at that time!");
                            }
                        }

                        footstep_definition.h_start = problem_data.velocity_constraint_problem_data.environment_surfaces->getSurfaceFromID(liftoff_surface_ID).surface_transform.translation()[2];
                        footstep_definition.h_end = problem_data.velocity_constraint_problem_data.environment_surfaces->getSurfaceFromID(touchdown_surface_ID).surface_transform.translation()[2];
                        footstep_definition.h_max = std::max(footstep_definition.h_start, footstep_definition.h_end) + problem_data.velocity_constraint_problem_data.ideal_offset_height;
                        footstep_definition.liftoff_time = liftoff_time;
                        footstep_definition.touchdown_time = touchdown_time;
                        double scaled = std::min(1.0, (touchdown_time - liftoff_time) / problem_data.velocity_constraint_problem_data.footstep_height_scaling);

                        casadi::SX mapped_t = (t - liftoff_time) / (touchdown_time - liftoff_time);
                        casadi::SX cfoot_vel = getFootstepVelocity(problem_data, ee.first);

                        casadi::SX ell_max_planar = problem_data.velocity_constraint_problem_data.max_following_leeway_planar;
                        casadi::SX ell_min_planar = problem_data.velocity_constraint_problem_data.min_following_leeway_planar;
                        casadi::SX ell_slope_planar = (ell_max_planar - ell_min_planar);
                        casadi::SX admissible_error_planar = ell_slope_planar * t + ell_min_planar;

                        casadi::Function lower_bound = casadi::Function("lower_bound", casadi::SXVector{t},
                                                                        casadi::SXVector{vertcat(casadi::SXVector{
                                                                            -admissible_error_planar,
                                                                            -admissible_error_planar})});
                        casadi::Function upper_bound = casadi::Function("upper_bound", casadi::SXVector{t},
                                                                        casadi::SXVector{vertcat(casadi::SXVector{
                                                                            admissible_error_planar,
                                                                            admissible_error_planar})});

                        G_vec.push_back(cfoot_vel(casadi::Slice(0, 2), 0));
                        lower_bound_vec.push_back(lower_bound(mapped_t).at(0));
                        upper_bound_vec.push_back(upper_bound(mapped_t).at(0));

                        casadi::Function footstep_function = casadi::Function("footstep_velocity", casadi::SXVector{t}, casadi::SXVector{createFootstepHeightFunction(t, footstep_definition, scaled)});
                        casadi::SX desired_velocity = footstep_function(casadi::SXVector{mapped_t}).at(0);

                        G_vec.push_back(cfoot_vel(2, 0));
                        lower_bound_vec.push_back(desired_velocity);
                        upper_bound_vec.push_back(desired_velocity);
                        if (ee.second->is_6d)
                        {
                            G_vec.push_back(cfoot_vel(casadi::Slice(3, 6), 0));
                            lower_bound_vec.push_back(casadi::SX::zeros(3, 1));
                            upper_bound_vec.push_back(casadi::SX::zeros(3, 1));
                        }
                    }
                    else
                    {
                        casadi::SX cfoot_vel = getFootstepVelocity(problem_data, ee.first);
                        if (ee.second->is_6d)
                        {
                            G_vec.push_back(cfoot_vel);
                            lower_bound_vec.push_back(casadi::SX::zeros(6, 1));
                            upper_bound_vec.push_back(casadi::SX::zeros(6, 1));
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
                constraint_data.lower_bound = casadi::Function("lower_bound", casadi::SXVector{t}, casadi::SXVector{vertcat(lower_bound_vec)});
                constraint_data.upper_bound = casadi::Function("upper_bound", casadi::SXVector{t}, casadi::SXVector{vertcat(upper_bound_vec)});

                constraint_data.metadata.name = "Velocity Constraint";
                int i = 0;
                for (auto ee : problem_data.velocity_constraint_problem_data.robot_end_effectors)
                {
                    if (mode[(*ee.second)])
                    {
                        // constraint_data.metadata.plot_titles.push_back("Stance Velocity Constraint " + ee.second->frame_name);
                        if (ee.second->is_6d)
                        {
                            constraint_data.metadata.plot_titles.push_back("Stance Velocity Constraint for " + ee.second->frame_name);
                            constraint_data.metadata.plot_groupings.push_back(std::make_tuple(i, i + 6));
                            constraint_data.metadata.plot_names.push_back({"Vx", "Vy", "Vz", "wx", "wy", "wz"});
                            i += 6;
                        }
                        else
                        {
                            constraint_data.metadata.plot_titles.push_back("Stance Velocity Constraint for " + ee.second->frame_name);
                            constraint_data.metadata.plot_groupings.push_back(std::make_tuple(i, i + 3));
                            constraint_data.metadata.plot_names.push_back({"Vx", "Vy", "Vz"});
                            i += 3;
                        }
                    }
                    else
                    {
                        constraint_data.metadata.plot_titles.push_back("Swing Velocity Constraint Vx for " + ee.second->frame_name);
                        constraint_data.metadata.plot_groupings.push_back(std::make_tuple(i, i + 1));
                        constraint_data.metadata.plot_names.push_back({"Vx"});
                        i += 1;
                        constraint_data.metadata.plot_titles.push_back("Swing Velocity Constraint Vy for " + ee.second->frame_name);
                        constraint_data.metadata.plot_groupings.push_back(std::make_tuple(i, i + 1));
                        constraint_data.metadata.plot_names.push_back({"Vy"});
                        i += 1;
                        constraint_data.metadata.plot_titles.push_back("Swing Velocity Constraint Vz for " + ee.second->frame_name);
                        constraint_data.metadata.plot_groupings.push_back(std::make_tuple(i, i + 1));
                        constraint_data.metadata.plot_names.push_back({"Vz"});
                        i += 1;
                    }
                }
            }
        }
    }
}