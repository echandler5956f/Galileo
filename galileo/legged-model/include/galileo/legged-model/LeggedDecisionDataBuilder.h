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
            struct LeggedDecisionProblemData
            {
                std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
                std::shared_ptr<contact::ContactSequence> contact_sequence;
                std::shared_ptr<legged::LeggedRobotStates> states;
                std::shared_ptr<opt::ADModel> ad_model;
                std::shared_ptr<opt::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;

                casadi::SX X0;
            };

            template <class ProblemData>
            class LeggedDecisionDataBuilder : public opt::DecisionDataBuilder<ProblemData>
            {
            public:
                LeggedDecisionDataBuilder() : opt::DecisionDataBuilder<ProblemData>() {}

                void buildDecisionData(const ProblemData &problem_data, int phase_index, opt::DecisionData &decision_data)
                {
                    contact::ContactMode mode = problem_data.legged_decision_problem_data.contact_sequence->getPhase(phase_index).mode;
                    std::shared_ptr<legged::LeggedRobotStates> states = problem_data.legged_decision_problem_data.states;

                    casadi::SX initial_guess_x = casadi::SX::zeros(states->nx, 1);
                    casadi::SX initial_guess_u = casadi::SX::zeros(states->nu, 1);

                    casadi::SX mass = problem_data.legged_decision_problem_data.ad_data->mass[0];
                    casadi::SX g = casadi::SX(9.81);
                    int num_in_contact = problem_data.legged_decision_problem_data.contact_sequence->numEndEffectorsInContactAtPhase(phase_index);

                    for (auto ee : problem_data.legged_decision_problem_data.robot_end_effectors)
                    {
                        if (mode[(*ee.second)])
                        {
                            initial_guess_u(std::get<0>(states->frame_id_to_index_range[ee.first]) + 2) = mass * g / num_in_contact;
                        }
                    }
                    initial_guess_x = problem_data.legged_decision_problem_data.X0;
                    decision_data.initial_guess = casadi::Function("DecisionInitialGuess", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{initial_guess_x, initial_guess_u});
                    decision_data.lower_bound = casadi::Function("DecisionLowerBounds", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{-casadi::inf * casadi::SX::ones(states->ndx, 1), -casadi::inf * casadi::SX::ones(states->nu, 1)});
                    decision_data.upper_bound = casadi::Function("DecisionUpperBounds", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{casadi::inf * casadi::SX::ones(states->ndx, 1), casadi::inf * casadi::SX::ones(states->nu, 1)});
                }
            };
        }
    }
}