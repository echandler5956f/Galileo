#pragma once

#include "galileo/opt/Solution.h"
#include "galileo/opt/Segment.h"
#include "galileo/opt/PseudospectralSegment.h"
#include "galileo/opt/PhaseSequence.h"
#include <chrono>

namespace galileo
{
    namespace opt
    {

        /**
         * @brief Callback function called at each iteration, used for debugging and plotting.
         * 
         */
        class IterationCallback : public casadi::Callback
        {
        public:
            /**
             * @brief Sparsity of the decision variables.
             * 
             */
            casadi::Sparsity x_;

            /**
             * @brief Sparsity of the objective function.
             * 
             */
            casadi::Sparsity f_;

            /**
             * @brief Sparsity of the constraints.
             * 
             */
            casadi::Sparsity g_;

            /**
             * @brief Sparsity of the Lagrange multipliers for the decision variables.
             * 
             */
            casadi::Sparsity lam_x_;

            /**
             * @brief Sparsity of the Lagrange multipliers for the constraints.
             * 
             */
            casadi::Sparsity lam_g_;

            /**
             * @brief Sparsity of the Lagrange multipliers for the parameters.
             * 
             */
            casadi::Sparsity lam_p_;

            /**
             * @brief Construct a new Iteration Callback object.
             * 
             */
            IterationCallback()
            {
                construct("f");
            }

            /**
             * @brief Destroy the Iteration Callback object.
             * 
             */
            ~IterationCallback() override {}

            /**
             * @brief Set the sparsity of the decision variables, objective function, constraints, Lagrange multipliers for the decision variables, constraints, and parameters.
             * 
             * @param x Sparsity of the decision variables
             * @param f Sparsity of the objective function
             * @param g Sparsity of the constraints
             * @param lam_x Sparsity of the Lagrange multipliers for the decision variables
             * @param lam_g Sparsity of the Lagrange multipliers for the constraints
             * @param lam_p Sparsity of the Lagrange multipliers for the parameters
             */
            void set_sparsity(casadi::Sparsity x, casadi::Sparsity f, casadi::Sparsity g, casadi::Sparsity lam_x, casadi::Sparsity lam_g, casadi::Sparsity lam_p)
            {
                this->x_ = x;
                this->f_ = f;
                this->g_ = g;
                this->lam_x_ = lam_x;
                this->lam_g_ = lam_g;
                this->lam_p_ = lam_p;
            }

            /**
             * @brief Initialize the callback function.
             * 
             */
            void init() override
            {
            }

            /**
             * @brief Get the number of inputs.
             * 
             * @return casadi_int The number of inputs
             */
            casadi_int get_n_in() override { return 6; }

            /**
             * @brief Get the number of outputs.
             * 
             * @return casadi_int The number of outputs
             */
            casadi_int get_n_out() override { return 1; }

            /**
             * @brief Evaluate the callback function.
             * 
             * @param arg The arguments to the callback function
             * @return std::vector<casadi::DM> The result of the callback function
             */
            std::vector<casadi::DM> eval(const std::vector<casadi::DM> &arg) const override
            {
                return std::vector<casadi::DM>{0};
            }

            /**
             * @brief Get the sparsity of the input at index i.
             * 
             * @param i The index of the input
             * @return casadi::Sparsity The sparsity of the input
             */
            casadi::Sparsity get_sparsity_in(casadi_int i) override
            {
                if (i == 0)
                    return x_;
                else if (i == 1)
                    return f_;
                else if (i == 2)
                    return g_;
                else if (i == 3)
                    return lam_x_;
                else if (i == 4)
                    return lam_g_;
                else if (i == 5)
                    return lam_p_;
                else
                    throw std::runtime_error("Invalid input index");
            }
        };

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
             * @param problem_ Problem data containing the objective function and dynamics
             * @param phase_sequence Sequence of phases including dynamics and timing information
             * @param builders_ Constraint builders used to build the constraints
             * @param decision_builder_ Decision builder used to build the decision data
             * @param opts_ Options to pass to the solver
             * @param nonlinear_solver_name_ Nonlinear solver name to use for the optimization
             */
            TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::shared_ptr<PhaseSequence<MODE_T>> phase_sequence, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders_, std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder_, casadi::Dict opts_, std::string nonlinear_solver_name_ = "ipopt");

            /**
             * @brief Initialize the finite elements.
             *
             * @param d The degree of the finite element polynomials
             * @param X0 The initial state to deviate from
             */
            void initFiniteElements(int d, casadi::DM X0);

            /**
             * @brief Optimize and return the solution.
             *
             * @return SXVector The solution
             */
            casadi::MXVector optimize();

            /**
             * @brief Collect the solution segments for each phase
             * 
             * @param times The times to evaluate the solution at
             * @param result Fill the result vector with the solution segments
             */
            void getSolutionSegments(const Eigen::VectorXd &times, std::vector<solution::solution_segment_data_t> &result);

            /**
             * @brief Get the constraint violations at the solution.
             * 
             * @param sol The solution to evaluate
             * @return std::vector<std::vector<solution::constraint_evaluations_t>> The constraint evaluations
             */
            std::vector<std::vector<solution::constraint_evaluations_t>> getConstraintViolations(solution::solution_t &sol) const;

            /**
             * @brief Get the indices of the segment times.
             *
             * @param times The times to evaluate the solution at
             * @param initial_time The initial time of the segment
             * @param end_time The end time of the segment
             * @return tuple_size_t The segment indices
             */
            tuple_size_t getSegmentIndices(Eigen::VectorXd times, double initial_time, double end_time) const;

        private:
            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<Segment>> trajectory;

            /**
             * @brief The ranges of the segment times.
             *
             */
            std::vector<tuple_size_t> segment_times_ranges;

            /**
             * @brief The solution vector.
             *
             */
            casadi::MXVector sol;

            /**
             * @brief Problem data containing constraints problem data.
             *
             */
            std::shared_ptr<ProblemData> problem;

            /**
             * @brief Problem data containing the objective function and dynamics.
             *
             */
            std::shared_ptr<GeneralProblemData> gp_data;

            /**
             * @brief Constraint builders used to build the constraints.
             *
             */
            std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders;

            /**
             * @brief Decision builder used to build the decision data.
             *
             */
            std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder;

            /**
             * @brief Constraint datas for each phase.
             *
             */
            std::vector<std::vector<std::shared_ptr<ConstraintData>>> constraint_datas_for_phase;

            /**
             * @brief Casadi solver options.
             *
             */
            casadi::Dict opts;

            /**
             * @brief Nonlinear solver to use.
             *
             */
            std::string nonlinear_solver_name;

            /**
             * @brief Nonlinear function solver.
             *
             */
            casadi::Function solver;

            /**
             * @brief Slicer to get the states.
             *
             */
            std::shared_ptr<States> state_indices;

            /**
             * @brief Sequence of phases including dynamics and timing information
             *
             */
            std::shared_ptr<PhaseSequence<MODE_T>> sequence;

            /**
             * @brief Range of decision variables per phase.
             *
             */
            std::vector<tuple_size_t> ranges_decision_variables;

            /**
             * @brief Vector of all decision variables.
             *
             */
            casadi::MXVector w;

            /**
             * @brief Vector of all constraint expressions.
             *
             */
            casadi::MXVector g;

            /**
             * @brief Vector of all general constraint lower bounds.
             *
             */
            std::vector<double> lbg;

            /**
             * @brief Vector of all general constraint upper bounds.
             *
             */
            std::vector<double> ubg;

            /**
             * @brief  Lower bounds associated with the decision variables.
             *
             */
            std::vector<double> lbw;

            /**
             * @brief Upper bounds associated with the decision variables.
             *
             */
            std::vector<double> ubw;

            /**
             * @brief Initial guess for the decision variables.
             */
            std::vector<double> w0;

            /**
             * @brief Expression for objective cost.
             *
             */
            casadi::MX J;

            /**
             * @brief Vector of all times where decision variables are evaluated.
             *
             */
            casadi::DM global_times;

            /**
             * @brief Callback function called at each iteration, used for debugging and plotting.
             *
             */
            IterationCallback callback;
        };

        template <class ProblemData, class MODE_T>
        TrajectoryOpt<ProblemData, MODE_T>::TrajectoryOpt(std::shared_ptr<ProblemData> problem_, std::shared_ptr<PhaseSequence<MODE_T>> phase_sequence, std::vector<std::shared_ptr<ConstraintBuilder<ProblemData>>> builders_, std::shared_ptr<DecisionDataBuilder<ProblemData>> decision_builder_, casadi::Dict opts_, std::string nonlinear_solver_name_)
        {
            this->problem = problem_;
            this->gp_data = problem_->gp_data;
            this->builders = builders_;
            this->decision_builder = decision_builder_;
            this->state_indices = problem_->states;
            this->sequence = phase_sequence;
            this->opts = opts_;
            this->nonlinear_solver_name = nonlinear_solver_name_;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::initFiniteElements(int d, casadi::DM X0)
        {
            assert(X0.size1() == state_indices->nx && X0.size2() == 1 && "Initial state must be a column vector");
            trajectory.clear();
            global_times = casadi::DM(0, 0);
            w.clear();
            g.clear();
            lbg.clear();
            ubg.clear();
            lbw.clear();
            ubw.clear();
            w0.clear();
            J = 0;

            casadi::MX prev_final_state = X0;
            casadi::MX prev_final_state_deviant;
            casadi::MX curr_initial_state_deviant;

            std::vector<double> equality_back_nx(state_indices->nx, 0.0);
            std::vector<double> equality_back_ndx(state_indices->ndx, 0.0);

            std::vector<std::shared_ptr<ConstraintData>> G;
            std::shared_ptr<DecisionData> Wdata = std::make_shared<DecisionData>();

            size_t num_phases = sequence->getNumPhases();

            printf("Starting initialization\n");
            auto begin_time = std::chrono::high_resolution_clock::now();
            auto start_time = std::chrono::high_resolution_clock::now();
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;

            for (size_t i = 0; i < num_phases; ++i)
            {
                start_time = std::chrono::high_resolution_clock::now();
                G.clear();
                for (auto builder : builders)
                {
                    std::shared_ptr<ConstraintData> con_data = std::make_shared<ConstraintData>();
                    builder->buildConstraint(*problem, i, *con_data);
                    G.push_back(con_data);
                }
                decision_builder->buildDecisionData(*problem, i, *Wdata);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for constraint building: " << elapsed.count() << std::endl;

                constraint_datas_for_phase.push_back(G);
                auto phase = sequence->phase_sequence_[i];

                start_time = std::chrono::high_resolution_clock::now();
                std::shared_ptr<Segment> segment = std::make_shared<PseudospectralSegment>(gp_data, phase.phase_dynamics, state_indices, d, phase.knot_points, phase.time_value / phase.knot_points);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for segment creation: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                segment->initializeSegmentTimeVector(global_times);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for local time vector initialization: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                segment->initializeInputTimeVector(global_times);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for input time vector initialization: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                segment->initializeKnotSegments(X0, prev_final_state);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for knot segment initialization: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                segment->initializeExpressionGraph(G, Wdata);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for expression graph initialization: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                segment->evaluateExpressionGraph(J, w, g);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for expression graph evaluation: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                ranges_decision_variables.push_back(segment->get_range_idx_decision_variables());

                trajectory.push_back(segment);

                segment->fill_lbg_ubg(lbg, ubg);
                segment->fill_lbw_ubw(lbw, ubw);
                segment->fill_w0(w0);
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for filling lbg, ubg, lbw, ubw, w0: " << elapsed.count() << std::endl;

                start_time = std::chrono::high_resolution_clock::now();
                /*Initial state constraint*/
                if (i == 0)
                {
                    auto curr_initial_state = segment->getInitialState();
                    g.push_back(prev_final_state - curr_initial_state);
                    lbg.insert(lbg.end(), equality_back_nx.begin(), equality_back_nx.end());
                    ubg.insert(ubg.end(), equality_back_nx.begin(), equality_back_nx.end());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    curr_initial_state_deviant = segment->getInitialStateDeviant();
                    /*For general jump map functions you can use the following syntax:*/
                    // g.push_back(jump_map_function(MXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    g.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    lbg.insert(lbg.end(), equality_back_ndx.begin(), equality_back_ndx.end());
                    ubg.insert(ubg.end(), equality_back_ndx.begin(), equality_back_ndx.end());
                }
                prev_final_state = segment->getFinalState();
                prev_final_state_deviant = segment->getFinalStateDeviant();

                /*Terminal cost*/
                if (i == num_phases - 1)
                {
                    J += gp_data->Phi(casadi::MXVector{prev_final_state}).at(0);
                }
                end_time = std::chrono::high_resolution_clock::now();
                elapsed = end_time - start_time;
                // std::cout << "Elapsed time for continuity and terminal cost: " << elapsed.count() << std::endl;
            }
            end_time = std::chrono::high_resolution_clock::now();
            elapsed = end_time - begin_time;
            std::cout << "Elapsed time for initialization: " << elapsed.count() << std::endl;
        }

        template <class ProblemData, class MODE_T>
        casadi::MXVector TrajectoryOpt<ProblemData, MODE_T>::optimize()
        {

            casadi::MXDict nlp = {{"x", vertcat(w)},
                                  {"f", J},
                                  {"g", vertcat(g)}};

            // callback.set_sparsity(vertcat(w).sparsity(), J.sparsity(), vertcat(g).sparsity(), vertcat(w).sparsity(), vertcat(g).sparsity(), casadi::Sparsity(0, 0));

            // opts["iteration_callback"] = callback;
            // opts["iteration_callback_step"] = 1; // Call the callback function at every iteration
            solver = casadi::nlpsol("solver", nonlinear_solver_name, nlp, opts);

            double time_from_funcs = 0.0;
            double time_just_solver = 0.0;
            casadi::DMDict arg;
            arg["lbg"] = lbg;
            arg["ubg"] = ubg;
            arg["lbx"] = lbw;
            arg["ubx"] = ubw;
            arg["x0"] = w0;
            casadi::DMDict result = solver(arg);
            w0 = result["x"].get_elements();
            casadi::Dict stats = solver.stats();
            if (nonlinear_solver_name == "ipopt")
            {
                time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_hess_l"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
                time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
                std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
                std::cout << "Total seconds from Ipopt w/o function: " << time_just_solver << std::endl;
            }
            else if (nonlinear_solver_name == "snopt")
            {
                time_from_funcs += (double)stats["t_wall_nlp_grad"] + (double)stats["t_wall_nlp_jac_f"] + (double)stats["t_wall_nlp_jac_g"];
                time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
                std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
                std::cout << "Total seconds from SNOPT w/o function: " << time_just_solver << std::endl;
            }
            auto full_sol = casadi::MX(result["x"]);

            for (size_t i = 0; i < trajectory.size(); ++i)
            {
                auto seg_sol = full_sol(casadi::Slice(0, casadi_int(std::get<1>(ranges_decision_variables[i]))));
                auto sol_i = trajectory[i]->extractSolution(seg_sol);
                if (i == 0)
                    sol = sol_i;
                else
                {
                    sol[0] = horzcat(sol[0], sol_i[0]);
                    sol[1] = horzcat(sol[1], sol_i[1]);
                }
            }
            return sol;
        }

        template <class ProblemData, class MODE_T>
        std::vector<std::vector<solution::constraint_evaluations_t>> TrajectoryOpt<ProblemData, MODE_T>::getConstraintViolations(solution::solution_t &result) const
        {
            auto clock_start_time = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<solution::constraint_evaluations_t>> constraint_evaluations;
            std::vector<solution::constraint_evaluations_t> phase_constraint_evaluations;

            casadi::DM dm_state_result = casadi::DM(casadi::Sparsity::dense(result.state_result.rows(), result.state_result.cols()));
            casadi::DM dm_input_result = casadi::DM(casadi::Sparsity::dense(result.input_result.rows(), result.input_result.cols()));
            casadi::DM dm_times = casadi::DM(result.times.rows(), result.times.cols());
            pinocchio::casadi::copy(result.state_result, dm_state_result);
            pinocchio::casadi::copy(result.input_result, dm_input_result);
            pinocchio::casadi::copy(result.times, dm_times);

            for (size_t i = 0; i < trajectory.size(); ++i)
            {
                phase_constraint_evaluations.clear();
                auto G = constraint_datas_for_phase[i];
                auto seg = trajectory[i];
                for (size_t j = 0; j < G.size(); ++j)
                {
                    std::shared_ptr<ConstraintData> con_data = G[j];
                    tuple_size_t seg_range = segment_times_ranges[i];
                    casadi_int start_idx = casadi_int(std::get<0>(seg_range));
                    casadi_int end_idx = casadi_int(std::get<1>(seg_range));

                    casadi::DM con_eval = con_data->G.map(end_idx - start_idx)(casadi::DMVector{
                                                                                   dm_state_result(casadi::Slice(0, dm_state_result.rows()), casadi::Slice(start_idx, end_idx)),
                                                                                   dm_input_result(casadi::Slice(0, dm_input_result.rows()), casadi::Slice(start_idx, end_idx))})
                                              .at(0);
                    casadi::DM con_lb = con_data->lower_bound.map(end_idx - start_idx)(casadi::DMVector{
                                                                                           dm_times(casadi::Slice(start_idx, end_idx), casadi::Slice(0, dm_times.columns()))})
                                            .at(0);
                    casadi::DM con_ub = con_data->upper_bound.map(end_idx - start_idx)(casadi::DMVector{
                                                                                           dm_times(casadi::Slice(start_idx, end_idx), casadi::Slice(0, dm_times.columns()))})
                                            .at(0);
                    Eigen::MatrixXd eval = Eigen::Map<Eigen::MatrixXd>(con_eval.get_elements().data(), con_eval.rows(), con_eval.columns()).transpose();
                    Eigen::MatrixXd lb = Eigen::Map<Eigen::MatrixXd>(con_lb.get_elements().data(), con_lb.rows(), con_lb.columns()).transpose();
                    Eigen::MatrixXd ub = Eigen::Map<Eigen::MatrixXd>(con_ub.get_elements().data(), con_ub.rows(), con_ub.columns()).transpose();

                    solution::constraint_evaluations_t con_evals;
                    con_evals.metadata = con_data->metadata;
                    con_evals.times = result.times.block(std::get<0>(seg_range), 0, std::get<1>(seg_range) - std::get<0>(seg_range), 1);
                    con_evals.evaluation = eval;
                    con_evals.lower_bounds = lb;
                    con_evals.upper_bounds = ub;
                    phase_constraint_evaluations.push_back(con_evals);
                }

                constraint_evaluations.push_back(phase_constraint_evaluations);
            }
            auto clock_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = clock_end_time - clock_start_time;
            std::cout << "Elapsed time to get constraint evaluations at solution: " << elapsed.count() << std::endl;
            return constraint_evaluations;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::getSolutionSegments(const Eigen::VectorXd &times, std::vector<solution::solution_segment_data_t> &result)
        {
            std::cout << "Getting solution segments" << std::endl;
            auto solx = sol[0];
            auto solu = sol[1];
            int state_count = 0;
            int input_count = 0;
            for (std::shared_ptr<Segment> seg : trajectory)
            {
                solution::solution_segment_data_t segment_data;

                std::shared_ptr<PseudospectralSegment> pseg = std::dynamic_pointer_cast<PseudospectralSegment>(seg);

                std::vector<double> state_times_vec = pseg->getSegmentTimes().get_elements();
                std::vector<double> input_times_vec = pseg->getUSegmentTimes().get_elements();
                segment_data.state_times = Eigen::Map<Eigen::VectorXd>(state_times_vec.data(), state_times_vec.size());
                segment_data.input_times = Eigen::Map<Eigen::VectorXd>(input_times_vec.data(), input_times_vec.size());
                segment_data.initial_time = state_times_vec[0];
                segment_data.end_time = state_times_vec[state_times_vec.size() - 1];

                segment_data.state_degree = pseg->getStateDegree();
                segment_data.input_degree = pseg->getInputDegree();
                segment_data.num_knots = pseg->getKnotNum();

                std::vector<double> solx_vec = casadi::MX::evalf(solx(casadi::Slice(0, solx.rows()), casadi::Slice(state_count, state_count + (pseg->getStateDegree() + 1) * pseg->getKnotNum()))).get_elements();
                std::vector<double> solu_vec = casadi::MX::evalf(solu(casadi::Slice(0, solu.rows()), casadi::Slice(input_count, input_count + (pseg->getInputDegree() + 1) * pseg->getKnotNum()))).get_elements();
                segment_data.solx_segment = Eigen::Map<Eigen::MatrixXd>(solx_vec.data(), solx.rows(), solx.columns());
                segment_data.solu_segment = Eigen::Map<Eigen::MatrixXd>(solu_vec.data(), solu.rows(), solu.columns());

                segment_data.state_poly = *(pseg->get_dXPoly());
                segment_data.input_poly = *(pseg->get_UPoly());

                state_count += (pseg->getStateDegree() + 1) * pseg->getKnotNum();
                input_count += (pseg->getInputDegree() + 1) * pseg->getKnotNum();

                tuple_size_t segment_indices = getSegmentIndices(times, segment_data.initial_time, segment_data.end_time);
                segment_times_ranges.push_back(segment_indices);

                result.push_back(segment_data);
            }
            std::cout << "Finished getting solution segments" << std::endl;
        }

        template <class ProblemData, class MODE_T>
        tuple_size_t TrajectoryOpt<ProblemData, MODE_T>::getSegmentIndices(Eigen::VectorXd times, double initial_time, double end_time) const
        {
            auto start_it = std::lower_bound(times.data(), times.data() + times.size(), initial_time);
            auto end_it = std::upper_bound(times.data(), times.data() + times.size(), end_time);
            return std::make_tuple(std::distance(times.data(), start_it), std::distance(times.data(), end_it));
        }
    }
}