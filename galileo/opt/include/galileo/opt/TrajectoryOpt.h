#pragma once

#include "galileo/opt/Solution.h"
#include "galileo/opt/PhaseSequence.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <map>
#include <set>
#include <algorithm>
#include <chrono>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief The trajectory optimization class. This class is responsible for
            initializing the finite elements, and optimizing the trajectory.
         *
         */
        template <class ProblemData, class MODE_T>
        class TrajectoryOpt
        {
        public:
            /**
             * @brief Construct a new Trajectory Opt object.
             *
             * @param problem Problem data containing the objective function and dynamics
             * @param phase_sequence Sequence of phases including dynamics and timing information
             * @param builders Constraint builders used to build the constraints
             * @param decision_builder Decision builder used to build the decision data
             * @param opts Options to pass to the solver
             * @param nonlinear_solver_name Nonlinear solver name to use for the optimization
             */
            TrajectoryOpt(std::shared_ptr<ProblemData> problem, std::shared_ptr<PhaseSequence<MODE_T>> phase_sequence, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders, std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder, casadi::Dict opts, std::string nonlinear_solver_name = "ipopt");

            /**
             * @brief Initialize the finite elements.
             *
             * @param d The degree of the finite element polynomials
             * @param horizon The time horizon for evaluating the finite elements.
             * If horizon < 0, the horizon is set to the total phase time.
             * This parameter determines which combinations of segments are converted into an nlp.
             * @param casadi_opts Casadi options
             */
            void InitFiniteElements(int d, double horizon = -1, casadi::Dict casadi_opts = casadi::Dict());

            /**
             * @brief Advance the finite elements.
             *
             * @param global_time The current global time
             * @param X0 The new initial state
             * @param Xf The new final state
             * @param w0 The new initial guess
             */
            void AdvanceFiniteElements(double global_time, casadi::DM X0, casadi::DM Xf, casadi::DM w0 = casadi::DM::zeros(0, 0));

            /**
             * @brief Optimize and return the solution.
             *
             * @return DMVector The solution
             */
            casadi::DMVector Optimize();

            /**
             * @brief Collect the solution segments for each phase
             *
             * @return std::vector<solution::solution_segment_data_t> The solution segments
             */
            std::vector<solution::solution_segment_data_t> getSolutionSegments() const;

            /**
             * @brief Get the constraint data segments for each phase.
             *
             * @return std::vector<std::vector<ConstraintData>> The constraint data segments
             */
            std::vector<std::vector<ConstraintData>> getConstraintDataSegments() const;

            /**
             * @brief Get the stacked decision variable times.
             *
             * @return casadi::DMVector The stacked decision variable times
             */
            casadi::DMVector getStackedDecisionVariableTimes() const;

        private:
            /**
             * @brief Reset the bounds and initial guess data.
             *
             */
            void ResetNLPInputData();

            /**
             * @brief Fill the ordered sequence ranges.
             *
             * For a horizon and a phase sequence, get the set of sequences of phases which could potentially occur during the moving horizon.
             * We will need an nlpsol object for each of these sequences.For instance, lets say we have a phase sequence of [0, 1, 2, 3, 4]
             * with phase periods of [0.25, 0.2, 0.25, 0.3, 0.25], and the user requests a horizon of 0.5.
             *
             * Then we need an nlpsol object for the following sequences of phases:
             *
             * [0, 1, 2], [1, 2, 3], [2, 3], [3, 4], [4]
             *
             * This function fills the ordered_sequence_ranges_ vector with the start and end indices of the sequences of phases.
             *
             */
            void FillOrderedSequenceRanges();

            /**
             * @brief Find the range of the sequence of phases which corresponds to the global time.
             *
             * @param global_time The global time
             * @param cumulative_times The cumulative times vector for the phases
             */
            tuple_size_t FindRangeFromTime(double global_time, const std::vector<double> &cumulative_times) const;

            /**
             * @brief Find the sequence range index which corresponds to the global time, given the moving horizon.
             *
             * @param global_time The global time
             * @return size_t The index of the sequence range
             */
            size_t FindSequenceRangeIndex(double global_time) const;

            /**
             * @brief Create the pseudospectral segments for the trajectory.
             *
             * @param X0_sym The initial state symbol
             * @param Xf_sym The final state symbol
             * @param d The degree of the pseudospectral segments
             * @param casadi_opts Casadi options
             *
             */
            void CreateSegments(casadi::MX X0_sym, casadi::MX Xf_sym, int d, casadi::Dict casadi_opts = casadi::Dict());

            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<PseudospectralSegment>> trajectory_;

            /**
             * @brief The solution vector.
             *
             */
            casadi::DMVector sol_;

            /**
             * @brief The previous solution.
             *
             */
            casadi::DM curr_solution_;

            /**
             * @brief Initial guess for the lagrange multipliers for bounds on X.
             *
             */
            casadi::DM lam_x0_;

            /**
             * @brief Initial guess for the lagrange multipliers for bounds on G.
             *
             */
            casadi::DM lam_g0_;

            /**
             * @brief Whether the solver is warm started.
             *
             */
            bool warm_ = false;

            /**
             * @brief Problem data containing constraints problem data.
             *
             */
            std::shared_ptr<ProblemData> problem_;

            /**
             * @brief Problem data containing the objective function and dynamics.
             *
             */
            std::shared_ptr<GeneralProblemData> gp_data_;

            /**
             * @brief Constraint builders used to build the constraints.
             *
             */
            std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders_;

            /**
             * @brief Decision builder used to build the decision data.
             *
             */
            std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder_;

            /**
             * @brief Constraint datas for each phase.
             *
             */
            std::vector<std::vector<ConstraintData>> constraint_datas_for_phase_;

            /**
             * @brief Casadi solver options.
             *
             */
            casadi::Dict opts_;

            /**
             * @brief Nonlinear solver to use.
             *
             */
            std::string nonlinear_solver_name_;

            /**
             * @brief Slicer to get the states.
             *
             */
            std::shared_ptr<States> state_indices_;

            /**
             * @brief Sequence of phases including dynamics and timing information
             *
             */
            std::shared_ptr<PhaseSequence<MODE_T>> sequence_;

            /**
             * @brief Range of the parameters shared between segments.
             *
             */
            tuple_size_t range_global_parameters_;

            /**
             * @brief Casadi likes to accumulate the function times, so we need to keep track of the previous time.
             *
             */
            double prev_time_from_funcs_ = 0.0;

            /**
             * @brief NLPProblemData struct containing decision variables, constraints, and objective function.
             *
             */
            NLPProblemData nlp_prob_data_;

            /**
             * @brief A vector of nonlinear programming problems.
             *
             */
            std::vector<NLP> nlps_;

            /**
             * @brief Index of the current nlp.
             *
             */
            size_t current_nlp_index_ = 0;

            /**
             * @brief Flag to indicate if the current_nlp_index_ has changed.
             *
             */
            bool nlp_changed_ = false;

            /**
             * @brief NLPInputData struct containing bounds and initial guess data.
             *
             */
            NLPInputData nlp_in_data_;

            /**
             * @brief For a horizon and a phase sequence, this is the set of sequences of segments which could potentially occur during the moving horizon
             *
             */
            std::vector<tuple_size_t> ordered_sequence_ranges_;

            /**
             * @brief The time horizon for evaluating the finite elements.
             * If horizon < 0, the horizon is set to the total phase time.
             * This parameter determines which combinations of segments are converted into an nlp.
             * Set from InitFiniteElements
             *
             */
            double horizon_ = -1;
        };

        // IMPORTANT: TrajectoryOpt assumes a fixed phase sequence.
        // This means we can handle periodic phase sequences, but not dynamic aperiodic phase sequences.
        // If you want dynamic aperiodic phase sequences, you need to create a new TrajectoryOpt object.
        template <class ProblemData, class MODE_T>
        TrajectoryOpt<ProblemData, MODE_T>::TrajectoryOpt(std::shared_ptr<ProblemData> problem, std::shared_ptr<PhaseSequence<MODE_T>> phase_sequence, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders, std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder, casadi::Dict opts, std::string nonlinear_solver_name)
        {
            problem_ = problem;
            gp_data_ = problem->gp_data;
            builders_ = builders;
            decision_builder_ = decision_builder;
            state_indices_ = problem->states;
            sequence_ = phase_sequence;
            opts_ = opts;
            nonlinear_solver_name_ = nonlinear_solver_name;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::InitFiniteElements(int d, double horizon, casadi::Dict casadi_opts)
        {
            std::cout << "Starting initialization" << std::endl;
            auto start = std::chrono::high_resolution_clock::now();

            trajectory_.clear();

            ResetNLPInputData();
            horizon_ = horizon;

            casadi::MX X0_sym = casadi::MX::sym("X0", state_indices_->nx, 1);
            casadi::MX Xf_sym = casadi::MX::sym("Xf", state_indices_->nx, 1);

            // All nlps use the same X0 and Xf symbolic parameters.
            range_global_parameters_ = tuple_size_t(0, 2 * state_indices_->nx);

            CreateSegments(X0_sym, Xf_sym, d, casadi_opts);
            FillOrderedSequenceRanges();

            nlps_.resize(ordered_sequence_ranges_.size());

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Creating pseudospectral segments took " << elapsed.count() * 1000. << " ms" << std::endl;

            start = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < ordered_sequence_ranges_.size(); ++i)
            {
                nlps_[i].nlp_prob_data.w.clear();
                nlps_[i].nlp_prob_data.g.clear();
                nlps_[i].nlp_prob_data.p.clear();
                nlps_[i].nlp_prob_data.J = 0;

                nlps_[i].nlp_prob_data.p.push_back(X0_sym);
                nlps_[i].nlp_prob_data.p.push_back(Xf_sym);

                casadi::MX prev_final_state_deviant;
                casadi::MXVector g_tmp_vector;

                tuple_size_t range = ordered_sequence_ranges_[i];

                // push back dummy vars so we can use the same indexing as the phases
                for (size_t j = 0; j < std::get<0>(range); ++j)
                {
                    nlps_[i].ranges_decision_variables.push_back(tuple_size_t(0, 0));
                    nlps_[i].ranges_parameters.push_back(tuple_size_t(0, 0));
                }

                for (size_t j = std::get<0>(range); j < std::get<1>(range); ++j)
                {
                    trajectory_[j]->FillNLPProblemData(nlps_[i].nlp_prob_data);

                    nlps_[i].ranges_decision_variables.push_back(trajectory_[j]->getRangeDecisionVariables());
                    nlps_[i].ranges_parameters.push_back(trajectory_[j]->getRangeParameters());

                    /*Initial state constraint*/
                    if (j == std::get<0>(range))
                    {
                        g_tmp_vector.push_back(X0_sym - trajectory_[j]->getInitialState());
                    }
                    /*Continuity constraint for the state deviant between phases*/
                    else if (j > std::get<0>(range))
                    {
                        g_tmp_vector.push_back(prev_final_state_deviant - trajectory_[j]->getInitialStateDeviant());
                    }
                    prev_final_state_deviant = trajectory_[j]->getFinalStateDeviant();

                    /*Terminal cost*/
                    if (j == std::get<1>(range) - 1)
                    {
                        nlps_[i].nlp_prob_data.J += gp_data_->Phi(casadi::MXVector{trajectory_[j]->getFinalState(), Xf_sym}).at(0);
                    }
                }
                nlps_[i].nlp_prob_data.g.push_back(vertcat(g_tmp_vector));

                casadi::MXDict nlp = {{"x", vertcat(nlps_[i].nlp_prob_data.w)},
                                      {"f", nlps_[i].nlp_prob_data.J},
                                      {"g", vertcat(nlps_[i].nlp_prob_data.g)},
                                      {"p", vertcat(nlps_[i].nlp_prob_data.p)}};

                nlps_[i].nlp_solver = casadi::nlpsol("solver" + std::to_string(i), nonlinear_solver_name_, nlp, opts_);
            }
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "Creating NLPs took " << elapsed.count() * 1000. << " ms" << std::endl;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::CreateSegments(casadi::MX X0_sym, casadi::MX Xf_sym, int d, casadi::Dict casadi_opts)
        {
            std::vector<ConstraintData> G;
            std::shared_ptr<DecisionData> Wdata = std::make_shared<DecisionData>();

            size_t num_phases = sequence_->getNumPhases();
            for (size_t i = 0; i < num_phases; ++i)
            {
                G.clear();
                for (auto builder : builders_)
                {
                    ConstraintData con_data;
                    builder->buildConstraint(*problem_, i, con_data);
                    G.push_back(con_data);
                }
                decision_builder_->buildDecisionData(*problem_, i, *Wdata);

                constraint_datas_for_phase_.push_back(G);
                auto phase = sequence_->getPhase(i);

                std::shared_ptr<PseudospectralSegment> segment = std::make_shared<PseudospectralSegment>(gp_data_, phase.phase_dynamics, phase.phase_cost, state_indices_, d, false, phase.knot_points);
                segment->InitializeExpressionGraph(G, Wdata, casadi_opts);
                segment->InitializeKnotSegments(X0_sym, Xf_sym);
                segment->EvaluateExpressionGraph();

                trajectory_.push_back(segment);
            }
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::FillOrderedSequenceRanges()
        {
            ordered_sequence_ranges_.clear();

            if (horizon_ > 0)
            {
                std::vector<double> cumulative_times = sequence_->getPhaseStartTimes();
                cumulative_times.push_back(cumulative_times.back() + sequence_->getPhase(sequence_->getNumPhases() - 1).time_value);

                // Define the interval for considering start times within each phase
                double interval = sequence_->getMinPeriod() / 7;
                for (double start_time = 0; start_time < cumulative_times.back(); start_time += interval)
                {
                    ordered_sequence_ranges_.push_back(FindRangeFromTime(start_time, cumulative_times));
                }

                ordered_sequence_ranges_.erase(std::unique(ordered_sequence_ranges_.begin(), ordered_sequence_ranges_.end()), ordered_sequence_ranges_.end());

                // Print the ranges
                for (tuple_size_t range : ordered_sequence_ranges_)
                {
                    std::cout << "(" << std::get<0>(range) << ", " << std::get<1>(range) << "), ";
                }
                std::cout << std::endl;
            }
            else
            {
                ordered_sequence_ranges_.push_back(tuple_size_t(0, sequence_->getNumPhases()));
            }
        }

        template <class ProblemData, class MODE_T>
        tuple_size_t TrajectoryOpt<ProblemData, MODE_T>::FindRangeFromTime(double global_time, const std::vector<double> &cumulative_times) const
        {
            if (global_time == cumulative_times.back())
            {
                return tuple_size_t(cumulative_times.size() - 2, cumulative_times.size() - 1);
            }
            auto start_it = std::upper_bound(cumulative_times.begin(), cumulative_times.end(), global_time);
            int start_index = std::distance(cumulative_times.begin(), start_it) - 1;
            if (start_index < 0)
            {
                std::cerr << "Error: start_index is less than 0" << std::endl;
            }

            double end_time = global_time + horizon_;
            auto end_it = std::upper_bound(cumulative_times.begin(), cumulative_times.end(), end_time);
            int end_index = std::distance(cumulative_times.begin(), end_it) - 1;
            end_index = std::min(end_index + 1, (int)cumulative_times.size() - 1);
            if (end_index < 0)
            {
                std::cerr << "Error: end_index is less than 0" << std::endl;
            }
            return tuple_size_t(start_index, end_index);
        }

        template <class ProblemData, class MODE_T>
        size_t TrajectoryOpt<ProblemData, MODE_T>::FindSequenceRangeIndex(double global_time) const
        {
            std::vector<double> cumulative_times = sequence_->getPhaseStartTimes();
            cumulative_times.push_back(cumulative_times.back() + sequence_->getPhase(sequence_->getNumPhases() - 1).time_value);

            tuple_size_t sequence_range = FindRangeFromTime(global_time, cumulative_times);

            // Use std::find to find the tuple in ordered_sequence_ranges_
            auto it = std::find(ordered_sequence_ranges_.begin(), ordered_sequence_ranges_.end(), sequence_range);

            // If the tuple was not found, return -1
            if (it == ordered_sequence_ranges_.end())
            {
                std::cerr << "Could not find the sequence range for the given global time." << std::endl;
                return -1;
            }

            // Otherwise, return the index of the found tuple
            return std::distance(ordered_sequence_ranges_.begin(), it);
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::AdvanceFiniteElements(double global_time, casadi::DM X0, casadi::DM Xf, casadi::DM w0)
        {
            auto start = std::chrono::high_resolution_clock::now();
            ResetNLPInputData();

            nlp_in_data_.p0.push_back(X0);
            nlp_in_data_.p0.push_back(Xf);

            casadi::DMVector lbg_tmp_vector, ubg_tmp_vector;

            size_t new_nlp_index;

            if (horizon_ < 0)
            {
                new_nlp_index = 0;
            }
            else
            {
                new_nlp_index = FindSequenceRangeIndex(global_time);
            }

            if (new_nlp_index != current_nlp_index_)
            {
                current_nlp_index_ = new_nlp_index;
                nlp_changed_ = true;
            }
            else
            {
                nlp_changed_ = false;
            }

            std::vector<double> cumulative_times = sequence_->getPhaseStartTimes();

            for (size_t i = std::get<0>(ordered_sequence_ranges_[current_nlp_index_]); i < std::get<1>(ordered_sequence_ranges_[current_nlp_index_]); ++i)
            {
                auto phase = sequence_->getPhase(i);
                double phase_start_time = sequence_->getPhaseStartTimes()[i];
                double phase_end_time = phase_start_time + phase.time_value;

                if (global_time > phase_end_time + 1e-6)
                {
                    std::cerr << "Global time is greater than the end of the phase. This should not be possible." << std::endl;
                    continue;
                }

                double update_time;
                double time_value;

                // Check if the global time is within the current phase
                if (phase_start_time < global_time && global_time < phase_start_time + phase.time_value)
                {
                    // Global time is within current phase
                    update_time = global_time;

                    // Calculate the remaining time in the phase
                    double remaining_time_in_phase = phase.time_value - (global_time - phase_start_time);

                    // Time value is the remaining time in the phase
                    time_value = remaining_time_in_phase;
                }
                else if (global_time + horizon_ > phase_start_time && global_time + horizon_ < phase_end_time && horizon_ > 0)
                {
                    // The current phase is ahead of the global time, but the phase is within the horizon
                    update_time = phase_start_time;

                    // The amount of time that is within the horizon
                    time_value = global_time + horizon_ - phase_start_time;
                }
                else
                {
                    // Global time is not within current phase
                    update_time = phase_start_time;

                    // Time value is the total time of the phase
                    time_value = phase.time_value;
                }

                // Determine if we need to update NLP data
                bool update_nlp_data = true;
                if (warm_ && w0.size1() > 0 && !nlp_changed_)
                {
                    update_nlp_data = false;
                }

                std::cout << "Phase: " << i << " covers " << time_value << " of the horizon at t = " << global_time << std::endl;

                // Update the trajectory and NLP data
                trajectory_[i]->Update(update_time, time_value / phase.knot_points, X0, Xf, update_nlp_data);
                trajectory_[i]->UpdateNLPInputData(nlp_in_data_, update_nlp_data);

                // If warm start is enabled and w0 has elements, add w0 to the NLP data, unless the nlp index has changed (in which case we are in a new sequence of phases)
                if (!update_nlp_data)
                {
                    nlp_in_data_.w0.push_back(w0);
                }

                /*Initial state constraint*/
                if (i == std::get<0>(ordered_sequence_ranges_[current_nlp_index_]))
                {
                    lbg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->nx, 1));
                    ubg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->nx, 1));
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > std::get<0>(ordered_sequence_ranges_[current_nlp_index_]))
                {
                    lbg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->ndx, 1));
                    ubg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->ndx, 1));
                }
            }
            nlp_in_data_.lbg.push_back(vertcat(lbg_tmp_vector));
            nlp_in_data_.ubg.push_back(vertcat(ubg_tmp_vector));

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Time to update finite elements: " << elapsed.count() * 1000. << " ms" << std::endl;

            Optimize();
        }

        template <class ProblemData, class MODE_T>
        casadi::DMVector TrajectoryOpt<ProblemData, MODE_T>::Optimize()
        {
            auto start = std::chrono::high_resolution_clock::now();
            casadi::DMDict arg;
            arg["lbg"] = vertcat(nlp_in_data_.lbg).get_elements();
            arg["ubg"] = vertcat(nlp_in_data_.ubg).get_elements();
            arg["lbx"] = vertcat(nlp_in_data_.lbw).get_elements();
            arg["ubx"] = vertcat(nlp_in_data_.ubw).get_elements();
            arg["p"] = vertcat(nlp_in_data_.p0).get_elements();
            arg["x0"] = vertcat(nlp_in_data_.w0).get_elements();

            if (warm_ && !nlp_changed_)
            {
                arg["lam_x0"] = lam_x0_;
                arg["lam_g0"] = lam_g0_;
            }
            else
            {
                warm_ = true;
            }

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Time to fill arguments: " << elapsed.count() * 1000. << " ms" << std::endl;

            start = std::chrono::high_resolution_clock::now();
            casadi::DMDict result = nlps_[current_nlp_index_].nlp_solver(arg);
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "Time to optimize: " << elapsed.count() * 1000. << " ms" << std::endl;

            start = std::chrono::high_resolution_clock::now();

            casadi::Dict stats = nlps_[current_nlp_index_].nlp_solver.stats();
            double time_from_funcs = 0.0, time_just_solver = 0.0;
            prev_time_from_funcs_ = nlp_changed_ ? 0.0 : prev_time_from_funcs_;
            if (nonlinear_solver_name_ == "ipopt")
            {
                time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"] - prev_time_from_funcs_;
                if (stats.find("t_wall_nlp_hess_l") != stats.end())
                    time_from_funcs += (double)stats["t_wall_nlp_hess_l"];
                time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
                std::cout << "Total seconds from Casadi functions: " << time_from_funcs * 1000. << " ms" << std::endl;
                std::cout << "Total seconds from Ipopt w/o function: " << time_just_solver * 1000. << " ms" << std::endl;
            }
            else if (nonlinear_solver_name_ == "snopt")
            {
                time_from_funcs += (double)stats["t_wall_nlp_grad"] + (double)stats["t_wall_nlp_jac_f"] + (double)stats["t_wall_nlp_jac_g"] - prev_time_from_funcs_;
                time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
                std::cout << "Total seconds from Casadi functions: " << time_from_funcs * 1000. << " ms" << std::endl;
                std::cout << "Total seconds from SNOPT w/o function: " << time_just_solver * 1000. << " ms" << std::endl;
            }
            prev_time_from_funcs_ += time_from_funcs;
            curr_solution_ = result["x"];
            lam_x0_ = result["lam_x"];
            lam_g0_ = result["lam_g"];
            casadi::DM p_sol = arg["p"];
            casadi::DM global_p_sol = p_sol(casadi::Slice(casadi_int(std::get<0>(range_global_parameters_)), casadi_int(std::get<1>(range_global_parameters_))));

            for (size_t i = std::get<0>(ordered_sequence_ranges_[current_nlp_index_]); i < std::get<1>(ordered_sequence_ranges_[current_nlp_index_]); ++i)
            {
                casadi::DM seg_sol = curr_solution_(casadi::Slice(casadi_int(std::get<0>(nlps_[current_nlp_index_].ranges_decision_variables[i])), casadi_int(std::get<1>(nlps_[current_nlp_index_].ranges_decision_variables[i]))));
                casadi::DM seg_p_sol = p_sol(casadi::Slice(casadi_int(std::get<0>(nlps_[current_nlp_index_].ranges_parameters[i])), casadi_int(std::get<1>(nlps_[current_nlp_index_].ranges_parameters[i]))));
                casadi::DM full_seg_p_sol = vertcat(casadi::DMVector{global_p_sol, seg_p_sol});
                casadi::DMVector sol_i = trajectory_[i]->ExtractSolution(seg_sol, full_seg_p_sol);
                if (i == std::get<0>(ordered_sequence_ranges_[current_nlp_index_]))
                    sol_ = sol_i;
                else
                {
                    sol_[0] = horzcat(sol_[0], sol_i[0]);
                    sol_[1] = horzcat(sol_[1], sol_i[1]);
                }
            }
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "Time to extract solution: " << elapsed.count() * 1000. << " ms" << std::endl;
            return sol_;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::ResetNLPInputData()
        {
            nlp_in_data_.w0.clear();
            nlp_in_data_.p0.clear();
            nlp_in_data_.lbw.clear();
            nlp_in_data_.ubw.clear();
            nlp_in_data_.lbg.clear();
            nlp_in_data_.ubg.clear();
        }

        template <class ProblemData, class MODE_T>
        casadi::DMVector TrajectoryOpt<ProblemData, MODE_T>::getStackedDecisionVariableTimes() const
        {
            casadi::DMVector result;
            for (size_t i = std::get<0>(ordered_sequence_ranges_[current_nlp_index_]); i < std::get<1>(ordered_sequence_ranges_[current_nlp_index_]); ++i)
            {
                casadi::DMVector stacked_segment_times = trajectory_[i]->getSegmentDecisionVariableTimes();
                result.insert(result.end(), stacked_segment_times.begin(), stacked_segment_times.end());
            }
            return result;
        }

        template <class ProblemData, class MODE_T>
        std::vector<std::vector<ConstraintData>> TrajectoryOpt<ProblemData, MODE_T>::getConstraintDataSegments() const
        {
            return constraint_datas_for_phase_;
        }

        template <class ProblemData, class MODE_T>
        std::vector<solution::solution_segment_data_t> TrajectoryOpt<ProblemData, MODE_T>::getSolutionSegments() const
        {
            std::vector<solution::solution_segment_data_t> result;
            auto solx = sol_[0];
            auto solu = sol_[1];
            int state_count = 0;
            int input_count = 0;
            for (size_t i = std::get<0>(ordered_sequence_ranges_[current_nlp_index_]); i < std::get<1>(ordered_sequence_ranges_[current_nlp_index_]); ++i)
            {
                std::shared_ptr<PseudospectralSegment> pseg = trajectory_[i];

                solution::solution_segment_data_t segment_data;

                std::vector<double> state_times_vec = pseg->getStateSegmentTimes().get_elements();
                std::vector<double> input_times_vec = pseg->getInputSegmentTimes().get_elements();
                segment_data.state_times = Eigen::Map<Eigen::VectorXd>(state_times_vec.data(), state_times_vec.size());
                segment_data.input_times = Eigen::Map<Eigen::VectorXd>(input_times_vec.data(), input_times_vec.size());
                segment_data.initial_time = state_times_vec[0];
                segment_data.end_time = state_times_vec[state_times_vec.size() - 1];

                segment_data.state_degree = pseg->getStateDegree();
                segment_data.input_degree = pseg->getInputDegree();
                segment_data.num_knots = pseg->getKnotNum();

                std::vector<double> solx_vec = solx(casadi::Slice(0, solx.rows()), casadi::Slice(state_count, state_count + (pseg->getStateDegree() + 1) * pseg->getKnotNum())).get_elements();
                std::vector<double> solu_vec = solu(casadi::Slice(0, solu.rows()), casadi::Slice(input_count, input_count + (pseg->getInputDegree() + 1) * pseg->getKnotNum())).get_elements();
                segment_data.solx_segment = Eigen::Map<Eigen::MatrixXd>(solx_vec.data(), solx.rows(), solx.columns());
                segment_data.solu_segment = Eigen::Map<Eigen::MatrixXd>(solu_vec.data(), solu.rows(), solu.columns());

                segment_data.state_poly = *(pseg->get_dXPoly());
                segment_data.input_poly = *(pseg->get_UPoly());

                state_count += (pseg->getStateDegree() + 1) * pseg->getKnotNum();
                input_count += (pseg->getInputDegree() + 1) * pseg->getKnotNum();

                result.push_back(segment_data);
            }
            return result;
        }
    }
}