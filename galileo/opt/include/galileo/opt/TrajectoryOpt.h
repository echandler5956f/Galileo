#pragma once

#include "galileo/opt/Solution.h"
#include "galileo/opt/PhaseSequence.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <map>
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
             */
            void InitFiniteElements(int d, double horizon = -1);

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
            std::vector<solution::solution_segment_data_t> getSolutionSegments();

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
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<PseudospectralSegment>> trajectory_;

            /**
             * @brief The ranges of the segment times.
             *
             */
            std::vector<tuple_size_t> segment_times_ranges_;

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
             * @brief Nonlinear function solver.
             *
             */
            casadi::Function solver_;

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
             * @brief Range of decision variables per phase.
             *
             */
            std::vector<tuple_size_t> ranges_decision_variables_;

            /**
             * @brief Range of parameters per phase.
             *
             */
            std::vector<tuple_size_t> ranges_segment_parameters_;

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
             * @brief NLPInputData struct containing bounds and initial guess data.
             * 
             */
            NLPInputData nlp_in_data_;

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
        void TrajectoryOpt<ProblemData, MODE_T>::InitFiniteElements(int d, double horizon)
        {
            trajectory_.clear();
            nlp_prob_data_.w.clear();
            nlp_prob_data_.g.clear();
            nlp_prob_data_.p.clear();
            nlp_prob_data_.J = 0;

            ResetNLPInputData();
            ranges_decision_variables_.clear();
            ranges_segment_parameters_.clear();

            horizon_ = horizon;

            casadi::MX prev_final_state_deviant;
            casadi::MXVector g_tmp_vector;

            std::vector<ConstraintData> G;
            std::shared_ptr<DecisionData> Wdata = std::make_shared<DecisionData>();

            size_t num_phases = sequence_->getNumPhases();
            std::cout << "Starting initialization" << std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            casadi::MX X0_sym = casadi::MX::sym("X0", state_indices_->nx, 1);
            casadi::MX Xf_sym = casadi::MX::sym("Xf", state_indices_->nx, 1);

            nlp_prob_data_.p.push_back(X0_sym);
            nlp_prob_data_.p.push_back(Xf_sym);
            range_global_parameters_ = tuple_size_t(0, 2 * state_indices_->nx);

            casadi::Dict casadi_opts = casadi::Dict();
            // {{"jit", true},
            //  {"jit_options.flags", "-Ofast -march=native -ffast-math"},
            //  {"jit_options.compiler", "gcc"},
            //  {"compiler", "shell"}};

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
                segment->FillNLPProblemData(nlp_prob_data_);

                ranges_decision_variables_.push_back(segment->getRangeDecisionVariables());
                ranges_segment_parameters_.push_back(segment->getRangeParameters());

                trajectory_.push_back(segment);

                /*Initial state constraint*/
                if (i == 0)
                {
                    g_tmp_vector.push_back(X0_sym - segment->getInitialState());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    /*For general jump map functions you can use the following syntax:*/
                    // g_tmp_vector.push_back(jump_map_function(MXVector{prev_final_state_deviant, segment->getInitialStateDeviant()}).at(0));
                    g_tmp_vector.push_back(prev_final_state_deviant - segment->getInitialStateDeviant());
                }
                prev_final_state_deviant = segment->getFinalStateDeviant();

                /*Terminal cost*/
                if (i == num_phases - 1)
                {
                    nlp_prob_data_.J += gp_data_->Phi(casadi::MXVector{segment->getFinalState(), Xf_sym}).at(0);
                }
            }
            // Ensures that the initial state constraint and continuity constraints are at the end of the constraint vector for convenience
            nlp_prob_data_.g.push_back(vertcat(g_tmp_vector));

            casadi::MXDict nlp = {{"x", vertcat(nlp_prob_data_.w)},
                                  {"f", nlp_prob_data_.J},
                                  {"g", vertcat(nlp_prob_data_.g)},
                                  {"p", vertcat(nlp_prob_data_.p)}};
            std::cout << "Finished initialization" << std::endl;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Time to initialize: " << elapsed.count() * 1000. << " ms" << std::endl;

            start = std::chrono::high_resolution_clock::now();
            solver_ = casadi::nlpsol("solver", nonlinear_solver_name_, nlp, opts_);
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "Time to create solver: " << elapsed.count() * 1000. << " ms" << std::endl;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::AdvanceFiniteElements(double global_time, casadi::DM X0, casadi::DM Xf, casadi::DM w0)
        {
            auto start = std::chrono::high_resolution_clock::now();
            ResetNLPInputData();

            nlp_in_data_.p0.push_back(X0);
            nlp_in_data_.p0.push_back(Xf);

            casadi::DMVector lbg_tmp_vector, ubg_tmp_vector;

            for (size_t i = 0; i < trajectory_.size(); ++i)
            {
                auto phase = sequence_->getPhase(i);
                if (global_time > sequence_->getPhaseStartTimes()[i] + phase.time_value + 1e-6)
                {
                    // handle removing this phase
                    std::cout << "Phase " << i << " has been removed" << std::endl;
                    continue;
                }

                // Calculate the start time of the current phase
                double phase_start_time = sequence_->getPhaseStartTimes()[i];

                // Check if the global time is within the current phase
                bool is_within_current_phase = (phase_start_time <= global_time) && (global_time < phase_start_time + phase.time_value);

                // Update time is either the global time (if within current phase) or the start time of the phase
                double update_time = is_within_current_phase ? global_time : phase_start_time;

                // Calculate the remaining time in the phase
                double remaining_time_in_phase = phase.time_value - (global_time - phase_start_time);

                // Time value depends on whether the update time is the global time
                double time_value = (update_time == global_time) ? remaining_time_in_phase : phase.time_value;

                // Determine if we need to update NLP data
                bool update_nlp_data = !(warm_ && w0.size1() > 0);

                // Update the trajectory and NLP data
                trajectory_[i]->Update(update_time, time_value / phase.knot_points, X0, Xf, update_nlp_data);
                trajectory_[i]->UpdateNLPInputData(nlp_in_data_, update_nlp_data);

                // If warm start is enabled and w0 has elements, add w0 to the NLP data
                if (warm_ && w0.size1() > 0)
                {
                    nlp_in_data_.w0.push_back(w0);
                }

                /*Initial state constraint*/
                if (i == 0)
                {
                    lbg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->nx, 1));
                    ubg_tmp_vector.push_back(casadi::DM::zeros(state_indices_->nx, 1));
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
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

            if (warm_)
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
            casadi::DMDict result = solver_(arg);
            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
            std::cout << "Time to optimize: " << elapsed.count() * 1000. << " ms" << std::endl;

            start = std::chrono::high_resolution_clock::now();
            casadi::Dict stats = solver_.stats();
            double time_from_funcs = 0.0;
            double time_just_solver = 0.0;
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

            for (size_t i = 0; i < trajectory_.size(); ++i)
            {
                casadi::DM seg_sol = curr_solution_(casadi::Slice(casadi_int(std::get<0>(ranges_decision_variables_[i])), casadi_int(std::get<1>(ranges_decision_variables_[i]))));
                casadi::DM seg_p_sol = p_sol(casadi::Slice(casadi_int(std::get<0>(ranges_segment_parameters_[i])), casadi_int(std::get<1>(ranges_segment_parameters_[i]))));
                casadi::DM full_seg_p_sol = vertcat(casadi::DMVector{global_p_sol, seg_p_sol});
                casadi::DMVector sol_i = trajectory_[i]->ExtractSolution(seg_sol, full_seg_p_sol);
                if (i == 0)
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
            for (size_t i = 0; i < trajectory_.size(); ++i)
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
        std::vector<solution::solution_segment_data_t> TrajectoryOpt<ProblemData, MODE_T>::getSolutionSegments()
        {
            std::vector<solution::solution_segment_data_t> result;
            auto solx = sol_[0];
            auto solu = sol_[1];
            int state_count = 0;
            int input_count = 0;
            for (std::shared_ptr<PseudospectralSegment> pseg : trajectory_)
            {
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