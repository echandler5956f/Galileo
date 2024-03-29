#pragma once

#include "galileo/legged-model/FrictionConeConstraintBuilder.h"
#include "galileo/legged-model/ContactConstraintBuilder.h"
#include "galileo/legged-model/VelocityConstraintBuilder.h"

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
                                       casadi::SX Xf)
                {
                    this->gp_data = gp_data_;
                    this->states = states_;

                    this->friction_cone_problem_data.environment_surfaces = environment_surfaces;
                    this->friction_cone_problem_data.contact_sequence = contact_sequence;
                    this->friction_cone_problem_data.states = states;
                    this->friction_cone_problem_data.ad_model = ad_model;
                    this->friction_cone_problem_data.ad_data = ad_data;
                    this->friction_cone_problem_data.robot_end_effectors = robot_end_effectors;
                    this->friction_cone_problem_data.x = x;
                    this->friction_cone_problem_data.u = u;
                    this->friction_cone_problem_data.t = t;
                    this->friction_cone_problem_data.mu = 0.7;
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

                    this->velocity_constraint_problem_data.ideal_offset_height = 0.08;
                    this->velocity_constraint_problem_data.max_following_leeway_planar = 0.25;
                    this->velocity_constraint_problem_data.min_following_leeway_planar = 1e-2;

                    this->velocity_constraint_problem_data.max_following_leeway_normal = 0.375;
                    this->velocity_constraint_problem_data.min_following_leeway_normal = this->velocity_constraint_problem_data.ideal_offset_height * 0.015;
                    this->velocity_constraint_problem_data.sigmoid_scaling = 50.;

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