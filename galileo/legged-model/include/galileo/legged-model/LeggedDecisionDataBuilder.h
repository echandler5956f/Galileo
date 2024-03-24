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

                // max height of CoM trajectory for initial guess
                double grf_reference = 100.0;
                double beta = -0.105;
                double gamma = 0.4;
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

                    // Gameplan: Generate joint and momentum trajectories using projectile motion (CoM trajectory).
                    // X-Y positions are found via linear interpolation between initial and final states.
                    // Z position is found using projectile motion of the nominal CoM offset from the base frame.
                    // Foot locations are estimated using the nominal joint configuration.

                    // If it is a flight phase, we use the ballistic CoM trajectory. If it is a normal phase, we use linear interpolation.

                    // std::vector<double> phase_timings = contact_sequence->getPhaseTiming();

                    // std::map<int, casadi::SX> initial_guess_x_per_phase;
                    // std::map<int, casadi::SX> initial_guess_u_per_phase;

                    // casadi::SX linear_interp = X0 + (Xf - X0) * problem_data.legged_decision_problem_data.t / t_f;
                    // casadi::SX quat1 = X0(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    // casadi::SX quat2 = Xf(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    // linear_interp(casadi::Slice(states->q_index + 3, states->q_index + states->nqb)) = math::quat_slerp(quat1, quat2, problem_data.legged_decision_problem_data.t / t_f);

                    // casadi::Function linear_interpolation = casadi::Function("LinearInterpolation", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{linear_interp});

                    // // // we are assuming we are starting from rest. this does not have to be the case in the future
                    // // casadi::SX q0_sym = casadi::SX::sym("q0", states->nq, 1);
                    // // ConfigVectorAD q0_vec = ConfigVectorAD::Zero(states->nq, 1);
                    // // pinocchio::casadi::copy(q0_sym, q0_vec);
                    // // ConfigVectorAD current_com_vec = pinocchio::centerOfMass(*(problem_data.legged_decision_problem_data.ad_model), *(problem_data.legged_decision_problem_data.ad_data), q0_vec, false);
                    // // casadi::SX current_com = casadi::SX::zeros(3, 1);
                    // // pinocchio::casadi::copy(current_com_vec, current_com);
                    // // casadi::SX com = casadi::Function("CoM", casadi::SXVector{q0_sym}, casadi::SXVector{current_com})(casadi::SXVector{casadi::SX(X0(casadi::Slice(states->q_index, states->q_index + states->nq)))}).at(0);

                    // casadi::SX nominal_com_offset = X0(casadi::Slice(states->q_index, states->q_index + 3));// - com;

                    // std::vector<int> non_flight_phase_indices;
                    // std::vector<int> flight_phase_indices;

                    // for (std::size_t i = 0; i < contact_sequence->phase_sequence_.size(); ++i)
                    // {
                    //     contact::ContactMode mode = contact_sequence->phase_sequence_[i].mode;
                    //     int num_in_contact = mode.numEndEffectorsInContact();
                    //     // if the first phase is a flight phase we treat it with linear interpolation
                    //     if (num_in_contact == 0 && i > 0)
                    //     {
                    //         flight_phase_indices.push_back(i);
                    //     }
                    //     else
                    //     {
                    //         non_flight_phase_indices.push_back(i);
                    //     }
                    // }

                    // // forward pass; linearly interpolate the states for the non-flight phases
                    // for (size_t i = 0; i < non_flight_phase_indices.size(); ++i)
                    // {
                    //     contact::ContactMode mode = contact_sequence->phase_sequence_[non_flight_phase_indices[i]].mode;
                    //     int num_in_contact = mode.numEndEffectorsInContact();
                    //     initial_guess_x_per_phase[non_flight_phase_indices[i]] = linear_interp;
                    //     casadi::SX initial_guess_u_phase = casadi::SX::zeros(states->nu, 1);

                    //     for (auto ee : problem_data.legged_decision_problem_data.robot_end_effectors)
                    //     {
                    //         if (mode[(*ee.second)])
                    //         {
                    //             initial_guess_u_phase(std::get<0>(states->frame_id_to_index_range[ee.first]) + 2) = mass * g / num_in_contact;
                    //         }
                    //     }
                    //     initial_guess_u_per_phase[non_flight_phase_indices[i]] = initial_guess_u_phase;
                    // }

                    // // backwards pass; get the non-flight phases that immediately precede the flight phases and generate the forces for a ballistic trajectory
                    // for (size_t i = 0; i < flight_phase_indices.size(); ++i)
                    // {
                    //     contact::ContactMode mode_before_flight = contact_sequence->phase_sequence_[flight_phase_indices[i] - 1].mode;
                    //     int num_in_contact_before_flight = mode_before_flight.numEndEffectorsInContact();
                    //     double t_start_before_flight = phase_timings[flight_phase_indices[i] - 1];
                    //     double t_lo = phase_timings[flight_phase_indices[i]];

                    //     casadi::SX t_scaled = (problem_data.legged_decision_problem_data.t - t_start_before_flight) / (t_lo - t_start_before_flight);

                    //     // this grf is just to get the vertical com acceleration, the actual guess for the contact forces comes later
                    //     casadi::SX grf_normal = problem_data.legged_decision_problem_data.grf_reference * num_in_contact_before_flight;
                    //     // this is the acceleration during the pre-flight phase
                    //     casadi::SX com_acceleration_normal = ((problem_data.legged_decision_problem_data.beta + (problem_data.legged_decision_problem_data.gamma * t_scaled)) * grf_normal / mass);

                    //     casadi::SX state_before_liftoff = linear_interpolation(casadi::SX{t_start_before_flight}).at(0);
                    //     casadi::SX com_position_normal_before_liftoff = state_before_liftoff(states->q_index + 2) - nominal_com_offset(2);
                    //     casadi::SX com_velocity_normal_before_liftoff = linear_interpolation(casadi::SX{t_start_before_flight}).at(0)(states->h_index + 2); // normalized centroidal momentum is simply the linear velocity of the com
                    //     casadi::SX com_position_normal_during_preflight = com_position_normal_before_liftoff + com_velocity_normal_before_liftoff * (problem_data.legged_decision_problem_data.t - t_start_before_flight) + 0.5 * com_acceleration_normal * (problem_data.legged_decision_problem_data.t - t_start_before_flight) * (problem_data.legged_decision_problem_data.t - t_start_before_flight);
                    //     casadi::SX base_position_z_during_preflight = com_position_normal_during_preflight + nominal_com_offset(2);

                    //     initial_guess_x_per_phase[flight_phase_indices[i] - 1] = linear_interp;
                    //     initial_guess_x_per_phase[flight_phase_indices[i] - 1](states->q_index + 2) = base_position_z_during_preflight;
                    //     initial_guess_u_per_phase[flight_phase_indices[i] - 1] = casadi::SX::zeros(states->nu, 1);
                    //     for (auto ee : problem_data.legged_decision_problem_data.robot_end_effectors)
                    //     {
                    //         if (mode_before_flight[(*ee.second)])
                    //         {
                    //             initial_guess_u_per_phase[flight_phase_indices[i] - 1](std::get<0>(states->frame_id_to_index_range[ee.first]) + 2) = (com_acceleration_normal + g) * mass / num_in_contact_before_flight;
                    //         }
                    //     }

                    //     // use the linear interpolation to get the com position up to the start of the flight phase, then integrate the acceleration to get the com velocity and position during the flight phase
                    //     casadi::SX state_at_liftoff = linear_interpolation(casadi::SX{t_lo}).at(0);
                    //     casadi::SX com_position_normal_at_liftoff = state_at_liftoff(states->q_index + 2) - nominal_com_offset(2);
                    //     casadi::SX com_velocity_normal_at_liftoff = linear_interpolation(casadi::SX{t_lo}).at(0)(states->h_index + 2); // normalized centroidal momentum is simply the linear velocity of the com
                    //     casadi::SX com_position_normal_during_flight = com_position_normal_at_liftoff + com_velocity_normal_at_liftoff * (problem_data.legged_decision_problem_data.t - t_lo) + 0.5 * com_acceleration_normal * (problem_data.legged_decision_problem_data.t - t_lo) * (problem_data.legged_decision_problem_data.t - t_lo);
                    //     casadi::SX base_position_z_during_flight = com_position_normal_during_flight + nominal_com_offset(2);

                    //     initial_guess_x_per_phase[flight_phase_indices[i]] = linear_interp;
                    //     initial_guess_x_per_phase[flight_phase_indices[i]](states->q_index + 2) = base_position_z_during_flight;
                    //     initial_guess_u_per_phase[flight_phase_indices[i]] = casadi::SX::zeros(states->nu, 1);
                    // }

                    // std::vector<casadi::SX> initial_guess_x_per_phase_vec;
                    // std::vector<casadi::SX> initial_guess_u_per_phase_vec;

                    // for (const auto &pair : initial_guess_x_per_phase)
                    // {
                    //     initial_guess_x_per_phase_vec.push_back(pair.second);
                    // }

                    // for (const auto &pair : initial_guess_u_per_phase)
                    // {
                    //     initial_guess_u_per_phase_vec.push_back(pair.second);
                    // }

                    // initial_guess_x = initial_guess_x_per_phase_vec[0];
                    // initial_guess_u = initial_guess_u_per_phase_vec[0];
                    // for (size_t i = 1; i < initial_guess_x_per_phase_vec.size(); i++)
                    // {
                    //     initial_guess_x = casadi::SX::if_else(problem_data.legged_decision_problem_data.t > phase_timings[i], initial_guess_x_per_phase_vec[i], initial_guess_x);
                    //     initial_guess_u = casadi::SX::if_else(problem_data.legged_decision_problem_data.t > phase_timings[i], initial_guess_u_per_phase_vec[i], initial_guess_u);
                    // }

                    // std::cout << "initial_guess_x: " << initial_guess_x << std::endl;
                    // std::cout << "initial_guess_u: " << initial_guess_u << std::endl;

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
                    initial_guess_x = X0 + (Xf -X0) * problem_data.legged_decision_problem_data.t / t_f;

                    casadi::SX quat1 = X0(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    casadi::SX quat2 = Xf(casadi::Slice(states->q_index + 3, states->q_index + states->nqb));
                    initial_guess_x(casadi::Slice(states->q_index + 3, states->q_index + states->nqb)) = math::quat_slerp(quat1, quat2, problem_data.legged_decision_problem_data.t / t_f);

                    decision_data.initial_guess = casadi::Function("DecisionInitialGuess", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{initial_guess_x, initial_guess_u});
                    decision_data.lower_bound = casadi::Function("DecisionLowerBounds", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{-casadi::inf * casadi::SX::ones(states->ndx, 1), -casadi::inf * casadi::SX::ones(states->nu, 1)});
                    decision_data.upper_bound = casadi::Function("DecisionUpperBounds", casadi::SXVector{problem_data.legged_decision_problem_data.t}, casadi::SXVector{casadi::inf * casadi::SX::ones(states->ndx, 1), casadi::inf * casadi::SX::ones(states->nu, 1)});
                }
            };
        }
    }
}