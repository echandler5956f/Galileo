#pragma once

#include "galileo/opt/Solution.h"
#include "galileo/opt/PhaseSequence.h"
#include "galileo/opt/PseudospectralSegment.h"
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
             * @brief Advance the finite elements.
             * 
             * @param time_advance The time to advance the finite elements
             * @param X0 The new initial state
             * @param phase The next phase (we check if we actually need the next phase based on the time_advance)
             */
            void advanceFiniteElements(double time_advance, casadi::DM X0, typename PhaseSequence<MODE_T>::Phase phase);

            /**
             * @brief Optimize and return the solution.
             *
             * @return SXVector The solution
             */
            casadi::MXVector optimize();

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

        private:
            /**
             * @brief Truncate part of or all of the first finite elements, such that the trajectory is advanced by time_advance.
             * 
             * @param time_advance The time to advance the finite elements
             * @param X0 The new initial state
             */
            void truncateFiniteElements(double time_advance, casadi::DM X0);

            /**
             * @brief A Trajectory is made up of segments of finite elements.
             *
             */
            std::vector<std::shared_ptr<PseudospectralSegment>> trajectory;

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
            std::vector<std::vector<ConstraintData>> constraint_datas_for_phase;

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
            casadi::DM local_times;
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
            local_times = casadi::DM(0, 0);
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
            casadi::MXVector g_tmp_vector;
            std::vector<double> lbg_tmp_vector;
            std::vector<double> ubg_tmp_vector;

            std::vector<ConstraintData> G;
            std::shared_ptr<DecisionData> Wdata = std::make_shared<DecisionData>();

            size_t num_phases = sequence->getNumPhases();
            std::cout << "Starting initialization" << std::endl;

            for (size_t i = 0; i < num_phases; ++i)
            {
                G.clear();
                for (auto builder : builders)
                {
                    ConstraintData con_data;
                    builder->buildConstraint(*problem, i, con_data);
                    G.push_back(con_data);
                }
                decision_builder->buildDecisionData(*problem, i, *Wdata);

                constraint_datas_for_phase.push_back(G);
                auto phase = sequence->getPhase(i);

                std::shared_ptr<PseudospectralSegment> segment = std::make_shared<PseudospectralSegment>(gp_data, phase.phase_dynamics, phase.phase_cost, state_indices, d, phase.knot_points, phase.time_value / phase.knot_points);
                segment->initializeSegmentTimeVector(local_times);
                segment->initializeInputTimeVector(local_times);
                segment->initializeKnotSegments(X0);
                segment->initializeExpressionGraph(G, Wdata);
                segment->evaluateExpressionGraph(J, w, g);

                ranges_decision_variables.push_back(segment->get_range_idx_decision_variables());

                segment->fill_lbg_ubg(lbg, ubg);
                segment->fill_lbw_ubw(lbw, ubw);
                segment->fill_w0(w0);

                trajectory.push_back(segment);

                /*Initial state constraint*/
                if (i == 0)
                {
                    auto curr_initial_state = segment->getInitialState();
                    g_tmp_vector.push_back(prev_final_state - curr_initial_state);
                    lbg_tmp_vector.insert(lbg_tmp_vector.end(), equality_back_nx.begin(), equality_back_nx.end());
                    ubg_tmp_vector.insert(ubg_tmp_vector.end(), equality_back_nx.begin(), equality_back_nx.end());
                }
                /*Continuity constraint for the state deviant between phases*/
                else if (i > 0)
                {
                    curr_initial_state_deviant = segment->getInitialStateDeviant();
                    /*For general jump map functions you can use the following syntax:*/
                    // g_tmp_vector.push_back(jump_map_function(MXVector{prev_final_state_deviant, curr_initial_state_deviant}).at(0));
                    g_tmp_vector.push_back(prev_final_state_deviant - curr_initial_state_deviant);
                    lbg_tmp_vector.insert(lbg_tmp_vector.end(), equality_back_ndx.begin(), equality_back_ndx.end());
                    ubg_tmp_vector.insert(ubg_tmp_vector.end(), equality_back_ndx.begin(), equality_back_ndx.end());
                }
                prev_final_state = segment->getFinalState();
                prev_final_state_deviant = segment->getFinalStateDeviant();

                /*Terminal cost*/
                if (i == num_phases - 1)
                {
                    J += gp_data->Phi(casadi::MXVector{prev_final_state}).at(0);
                }
            }
            // Ensures that the initial state constraint and continuity constraints are at the end of the constraint vector for convenience
            g.insert(g.begin(), g_tmp_vector.begin(), g_tmp_vector.end());
            lbg.insert(lbg.begin(), lbg_tmp_vector.begin(), lbg_tmp_vector.end());
            ubg.insert(ubg.begin(), ubg_tmp_vector.begin(), ubg_tmp_vector.end());
            std::cout << "Finished initialization" << std::endl;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::advanceFiniteElements(double time_advance, casadi::DM X0, typename PhaseSequence<MODE_T>::Phase phase)
        {
            truncateFiniteElements(time_advance, X0);
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::truncateFiniteElements(double time_advance, casadi::DM X0)
        {
            // For now we assume that we are advancing 1 phase, and that all phases have the same duration/knot number/degree. 
            // This is a big assumption that we will fix later.

            lbw.erase(lbw.begin() + std::get<0>(trajectory[0]->get_range_idx_decision_bounds()), lbw.begin() + std::get<1>(trajectory[0]->get_range_idx_decision_bounds()));
            ubw.erase( ubw.begin() + std::get<0>(trajectory[0]->get_range_idx_decision_bounds()), ubw.begin() + std::get<1>(trajectory[0]->get_range_idx_decision_bounds()));
            w0.erase(w0.begin() + std::get<0>(trajectory[0]->get_range_idx_decision_bounds()), w0.begin() + std::get<1>(trajectory[0]->get_range_idx_decision_bounds()));
            w.erase(w.begin() + std::get<0>(trajectory[0]->get_range_idx_decision_variables()), w.begin() + std::get<1>(trajectory[0]->get_range_idx_decision_variables()));

            g.erase(g.begin() + std::get<0>(trajectory[0]->get_range_idx_constraint_expressions()), g.begin() + std::get<1>(trajectory[0]->get_range_idx_constraint_expressions()));
            lbg.erase(lbg.begin() + std::get<0>(trajectory[0]->get_range_idx_constraint_bounds()), lbg.begin() + std::get<1>(trajectory[0]->get_range_idx_constraint_bounds()));
            ubg.erase(ubg.begin() + std::get<0>(trajectory[0]->get_range_idx_constraint_bounds()), ubg.begin() + std::get<1>(trajectory[0]->get_range_idx_constraint_bounds()));

            local_times = local_times(casadi::Slice(std::get<1>(trajectory[0]->get_range_idx_time()), local_times.size1()));
            local_times = local_times - local_times(0);

            trajectory.erase(trajectory.begin());

            g.push_back(X0 - trajectory[0]->getInitialStateDeviant());
            lbg.insert(lbg.end(), lbg.begin(), lbg.begin() + state_indices->ndx);
            ubg.insert(ubg.end(), ubg.begin(), ubg.begin() + state_indices->ndx);
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
                time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
                if (stats.find("t_wall_nlp_hess_l") != stats.end())
                    time_from_funcs += (double)stats["t_wall_nlp_hess_l"];
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
        std::vector<std::vector<ConstraintData>> TrajectoryOpt<ProblemData, MODE_T>::getConstraintDataSegments() const
        {
            return constraint_datas_for_phase;
        }

        template <class ProblemData, class MODE_T>
        std::vector<solution::solution_segment_data_t> TrajectoryOpt<ProblemData, MODE_T>::getSolutionSegments()
        {
            std::vector<solution::solution_segment_data_t> result;
            auto solx = sol[0];
            auto solu = sol[1];
            int state_count = 0;
            int input_count = 0;
            for (std::shared_ptr<PseudospectralSegment> pseg : trajectory)
            {
                solution::solution_segment_data_t segment_data;

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

                result.push_back(segment_data);
            }
            return result;
        }
    }
}