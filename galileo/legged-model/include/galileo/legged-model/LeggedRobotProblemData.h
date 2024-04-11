#pragma once

#include "galileo/legged-model/FrictionConeConstraintBuilder.h"
#include "galileo/legged-model/ContactConstraintBuilder.h"
#include "galileo/legged-model/NormalVelocityEqualityConstraintBuilder.h"
#include "galileo/legged-model/LeggedDecisionDataBuilder.h"

namespace galileo
{
    namespace legged
    {
        namespace constraints
        {
            class LeggedRobotProblemData
            {
            public:
                LeggedRobotProblemData(std::shared_ptr<opt::GeneralProblemData> gp_data_,
                                       std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces,
                                       std::shared_ptr<contact::ContactSequence> contact_sequence,
                                       std::shared_ptr<legged::LeggedRobotStates> states_,
                                       std::shared_ptr<legged::ADModel> ad_model,
                                       std::shared_ptr<legged::ADData> ad_data,
                                       contact::RobotEndEffectors robot_end_effectors,
                                       casadi::SX x,
                                       casadi::SX u,
                                       casadi::SX t,
                                       casadi::SX X0,
                                       casadi::SX Xf,
                                       JointLimits joint_limits,
                                       std::map<std::string, double> opts)
                {
                    this->gp_data = gp_data_;
                    this->states = states_;

                    phase_sequence = contact_sequence;

                    this->friction_cone_problem_data.environment_surfaces = environment_surfaces;
                    this->friction_cone_problem_data.contact_sequence = contact_sequence;
                    this->friction_cone_problem_data.states = states;
                    this->friction_cone_problem_data.ad_model = ad_model;
                    this->friction_cone_problem_data.ad_data = ad_data;
                    this->friction_cone_problem_data.robot_end_effectors = robot_end_effectors;
                    this->friction_cone_problem_data.x = x;
                    this->friction_cone_problem_data.u = u;
                    this->friction_cone_problem_data.t = t;
                    if (opts.find("mu") != opts.end())
                        this->friction_cone_problem_data.mu = opts["mu"];
                    if (opts.find("normal_force_max") != opts.end())
                        this->friction_cone_problem_data.normal_force_max = opts["normal_force_max"];
                    this->friction_cone_problem_data.approximation_order = FrictionConeProblemData::ApproximationOrder::FIRST_ORDER;

                    this->contact_constraint_problem_data.environment_surfaces = environment_surfaces;
                    this->contact_constraint_problem_data.contact_sequence = contact_sequence;
                    this->contact_constraint_problem_data.states = states;
                    this->contact_constraint_problem_data.ad_model = ad_model;
                    this->contact_constraint_problem_data.ad_data = ad_data;
                    this->contact_constraint_problem_data.robot_end_effectors = robot_end_effectors;
                    this->contact_constraint_problem_data.x = x;
                    this->contact_constraint_problem_data.u = u;
                    this->contact_constraint_problem_data.t = t;

                    this->velocity_constraint_problem_data.environment_surfaces = environment_surfaces;
                    this->velocity_constraint_problem_data.contact_sequence = contact_sequence;
                    this->velocity_constraint_problem_data.states = states;
                    this->velocity_constraint_problem_data.ad_model = ad_model;
                    this->velocity_constraint_problem_data.ad_data = ad_data;
                    this->velocity_constraint_problem_data.robot_end_effectors = robot_end_effectors;
                    this->velocity_constraint_problem_data.x = x;
                    this->velocity_constraint_problem_data.u = u;
                    this->velocity_constraint_problem_data.t = t;
                    if (opts.find("ideal_offset_height") != opts.end())
                        this->velocity_constraint_problem_data.ideal_offset_height = opts["ideal_offset_height"];
                    if (opts.find("footstep_height_scaling") != opts.end())
                        this->velocity_constraint_problem_data.footstep_height_scaling = opts["footstep_height_scaling"];
                    if (opts.find("max_following_leeway_planar") != opts.end())
                        this->velocity_constraint_problem_data.max_following_leeway_planar = opts["max_following_leeway_planar"];
                    if (opts.find("min_following_leeway_planar") != opts.end())
                        this->velocity_constraint_problem_data.min_following_leeway_planar = opts["min_following_leeway_planar"];
                    if (opts.find("footstep_vel_start") != opts.end())
                        this->velocity_constraint_problem_data.footstep_vel_start = opts["footstep_vel_start"];
                    if (opts.find("footstep_vel_end") != opts.end())
                        this->velocity_constraint_problem_data.footstep_vel_end = opts["footstep_vel_end"];

                    this->legged_decision_problem_data.environment_surfaces = environment_surfaces;
                    this->legged_decision_problem_data.contact_sequence = contact_sequence;
                    this->legged_decision_problem_data.states = states;
                    this->legged_decision_problem_data.ad_model = ad_model;
                    this->legged_decision_problem_data.ad_data = ad_data;
                    this->legged_decision_problem_data.robot_end_effectors = robot_end_effectors;
                    this->legged_decision_problem_data.x = x;
                    this->legged_decision_problem_data.u = u;
                    this->legged_decision_problem_data.t = t;
                    this->legged_decision_problem_data.X0 = X0;
                    this->legged_decision_problem_data.Xf = Xf;
                    this->legged_decision_problem_data.joint_limits = joint_limits;
                }
                std::shared_ptr<opt::PhaseSequence<contact::ContactMode>> phase_sequence;
                std::shared_ptr<opt::GeneralProblemData> gp_data;
                std::shared_ptr<legged::LeggedRobotStates> states;
                FrictionConeProblemData friction_cone_problem_data;
                ContactConstraintProblemData contact_constraint_problem_data;
                VelocityConstraintProblemData velocity_constraint_problem_data;

                LeggedDecisionProblemData legged_decision_problem_data;
            };
        }
    }
}