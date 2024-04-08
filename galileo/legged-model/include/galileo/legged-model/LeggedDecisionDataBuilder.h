#pragma once

#include "galileo/opt/Constraint.h"
#include "galileo/math/LieAlgebra.h"
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
                std::shared_ptr<legged::ADModel> ad_model;
                std::shared_ptr<legged::ADData> ad_data;
                contact::RobotEndEffectors robot_end_effectors;
                casadi::SX x;
                casadi::SX u;
                casadi::SX t;

                casadi::SX X0;
                casadi::SX Xf;
            };

            template <class ProblemData>
            class LeggedDecisionDataBuilder : public opt::DecisionDataBuilder<ProblemData>
            {
            public:
                LeggedDecisionDataBuilder() : opt::DecisionDataBuilder<ProblemData>() {}

                void buildDecisionData(const ProblemData &problem_data, int phase_index, opt::DecisionData &decision_data)
                {
                    std::shared_ptr<legged::LeggedRobotStates> states = problem_data.legged_decision_problem_data.states;
                    std::shared_ptr<contact::ContactSequence> contact_sequence = problem_data.legged_decision_problem_data.contact_sequence;

                    casadi::SX X0 = problem_data.legged_decision_problem_data.X0;
                    casadi::SX Xf = problem_data.legged_decision_problem_data.Xf;

                    casadi::SX t = problem_data.legged_decision_problem_data.t;
                    double t_f = problem_data.legged_decision_problem_data.contact_sequence->getDT();

                    casadi::SX initial_guess_x = casadi::SX::zeros(states->nx, 1);
                    casadi::SX initial_guess_u = casadi::SX::zeros(states->nu, 1);

                    casadi::SX mass = problem_data.legged_decision_problem_data.ad_data->mass[0];
                    casadi::SX g = casadi::SX(9.81);

                    for (auto ee : problem_data.legged_decision_problem_data.robot_end_effectors)
                    {
                        contact::ContactMode mode = contact_sequence->phase_sequence_[phase_index].mode;
                        int num_in_contact = mode.numEndEffectorsInContact();
                        if (mode[(*ee.second)])
                        {
                            initial_guess_u(std::get<0>(states->frame_id_to_index_range[ee.first]) + 2) = mass * g / num_in_contact;
                        }
                    }
                    // linearly interpolate between initial and final states
                    initial_guess_x = X0 + (Xf - X0) * t / t_f;

                    // // Calculate the average linear momentum required to reach the final state and use it as the initial guess
                    // casadi::SX q0 = states->get_q(X0);
                    // casadi::SX qf = states->get_q(Xf);
                    // casadi::SX avg_lin_velocity = (qf(casadi::Slice(0, 3)) - q0(casadi::Slice(0, 3))) / t_f; // the average linear velocity required to reach the final state
                    // // casadi::SX lin_velocity = avg_lin_velocity + casadi::SX::sin(2 * M_PI * t / t_f - M_PI / 2);
                    // initial_guess_x(casadi::Slice(states->h_index, states->h_index + 3)) = avg_lin_velocity;

                    // // apply this linear velocity to the guess for the base frame position
                    // // initial_guess_x(casadi::Slice(states->q_index, states->q_index + 3)) = q0(casadi::Slice(0, 3)) + avg_lin_velocity * t - avg_lin_velocity * (t_f / (2 * M_PI)) * casadi::SX::sin(2 * M_PI * t / t_f);

                    casadi::SX quat1 = X0(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    casadi::SX quat2 = Xf(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    initial_guess_x(casadi::Slice(states->q_index + 3, states->q_index + states->nqb)) = math::quat_slerp(quat1, quat2, t / t_f);

                    decision_data.initial_guess = casadi::Function("DecisionInitialGuess", casadi::SXVector{t}, casadi::SXVector{initial_guess_x, initial_guess_u});
                    decision_data.lower_bound = casadi::Function("DecisionLowerBounds", casadi::SXVector{t}, casadi::SXVector{-casadi::inf * casadi::SX::ones(states->ndx, 1), -casadi::inf * casadi::SX::ones(states->nu, 1)});
                    decision_data.upper_bound = casadi::Function("DecisionUpperBounds", casadi::SXVector{t}, casadi::SXVector{casadi::inf * casadi::SX::ones(states->ndx, 1), casadi::inf * casadi::SX::ones(states->nu, 1)});

                    // // print initial guess for 0 to t_f in 100 steps
                    // std::cout << "--------------------------------" << std::endl;
                    // for (int i = 0; i < 100; i++)
                    // {
                    //     double t_val = t_f * i / 100;
                    //     casadi::DMVector initial_guess = decision_data.initial_guess(casadi::DM(t_val));
                    //     std::cout << "t = " << t_val << std::endl;
                    //     casadi::DM x = initial_guess[0];
                    //     std::cout << "pos = " << x(casadi::Slice(states->q_index, states->q_index + 3)) << std::endl;
                    //     std::cout << "--------------------------------" << std::endl;
                    //     // std::cout << "u = " << initial_guess[1] << std::endl;
                    // }
                }
            };
        }
    }
}