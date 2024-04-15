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
             */
            void InitFiniteElements(int d);

            /**
             * @brief Advance the finite elements.
             *
             * @param X0 The new initial state
             * @param Xf The new final state
             */
            void AdvanceFiniteElements(casadi::DM X0, casadi::DM Xf);

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

        private:
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
            casadi::DM prev_solution_;

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
            std::vector<tuple_size_t> ranges_parameters_;

            /**
             * @brief NLPData struct containing decision variables, constraints, and objective function.
             *
             */
            NLPData nlp_data_;
        };

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
        void TrajectoryOpt<ProblemData, MODE_T>::InitFiniteElements(int d)
        {
            assert(X0.size1() == state_indices_->nx && X0.size2() == 1 && "Initial state must be a column vector");
            trajectory_.clear();
            nlp_data_.w.clear();
            nlp_data_.g.clear();
            nlp_data_.p.clear();
            nlp_data_.J = 0;

            nlp_data_.w0.clear();
            nlp_data_.p0.clear();
            nlp_data_.lbw.clear();
            nlp_data_.ubw.clear();
            nlp_data_.lbg.clear();
            nlp_data_.ubg.clear();
            ranges_decision_variables_.clear();
            ranges_parameters_.clear();

            casadi::MX prev_final_state_deviant;

            casadi::MXVector g_tmp_vector;

            std::vector<ConstraintData> G;
            std::shared_ptr<DecisionData> Wdata = std::make_shared<DecisionData>();

            size_t num_phases = sequence_->getNumPhases();
            std::cout << "Starting initialization" << std::endl;
            casadi::MX X0_sym = casadi::MX::sym("X0", state_indices_->nx, 1);
            casadi::MX Xf_sym = casadi::MX::sym("Xf", state_indices_->nx, 1);

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
                segment->InitializeExpressionGraph(G, Wdata);
                segment->InitializeKnotSegments(X0_sym, Xf_sym);
                segment->EvaluateExpressionGraph();
                segment->FillNLPData(nlp_data_);

                ranges_decision_variables_.push_back(segment->getRangeDecisionVariables());
                ranges_parameters_.push_back(segment->getRangeParameters());

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
                    nlp_data_.J += gp_data_->Phi(casadi::MXVector{segment->getFinalState(), Xf_sym}).at(0);
                }
            }
            // Ensures that the initial state constraint and continuity constraints are at the end of the constraint vector for convenience
            nlp_data_.g.push_back(vertcat(g_tmp_vector));

            casadi::MXDict nlp = {{"x", vertcat(nlp_data_.w)},
                                  {"f", nlp_data_.J},
                                  {"g", vertcat(nlp_data_.g)},
                                  {"p", vertcat(nlp_data_.p)}};

            solver_ = casadi::nlpsol("solver", nonlinear_solver_name_, nlp, opts_);
            std::cout << "Finished initialization" << std::endl;
        }

        template <class ProblemData, class MODE_T>
        void TrajectoryOpt<ProblemData, MODE_T>::AdvanceFiniteElements(casadi::DM X0, casadi::DM Xf)
        {
            nlp_data_.w0.clear();
            nlp_data_.p0.clear();
            nlp_data_.lbw.clear();
            nlp_data_.ubw.clear();
            nlp_data_.lbg.clear();
            nlp_data_.ubg.clear();

            casadi::DMVector lbg_tmp_vector, ubg_tmp_vector;

            for (size_t i = 0; i < trajectory_.size(); ++i)
            {
                auto phase = sequence_->getPhase(i);
                trajectory_[i]->Update(sequence_->getPhaseStartTimes()[i], phase.time_value / phase.knot_points, X0, Xf);
                trajectory_[i]->UpdateNLPData(nlp_data_);

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
            nlp_data_.lbg.push_back(vertcat(lbg_tmp_vector));
            nlp_data_.ubg.push_back(vertcat(ubg_tmp_vector));

            Optimize();
        }

        template <class ProblemData, class MODE_T>
        casadi::DMVector TrajectoryOpt<ProblemData, MODE_T>::Optimize()
        {
            // double time_from_funcs = 0.0;
            // double time_just_solver = 0.0;
            casadi::DMDict arg;
            arg["lbg"] = vertcat(nlp_data_.lbg).get_elements();
            arg["ubg"] = vertcat(nlp_data_.ubg).get_elements();
            arg["lbx"] = vertcat(nlp_data_.lbw).get_elements();
            arg["ubx"] = vertcat(nlp_data_.ubw).get_elements();
            arg["x0"] = vertcat(nlp_data_.w0).get_elements();
            arg["p"] = vertcat(nlp_data_.p0).get_elements();

            casadi::DMDict result = solver_(arg);
            // casadi::Dict stats = solver_.stats();
            // if (nonlinear_solver_name_ == "ipopt")
            // {
            //     time_from_funcs += (double)stats["t_wall_nlp_jac_g"] + (double)stats["t_wall_nlp_grad_f"] + (double)stats["t_wall_nlp_g"] + (double)stats["t_wall_nlp_f"];
            //     if (stats.find("t_wall_nlp_hess_l") != stats.end())
            //         time_from_funcs += (double)stats["t_wall_nlp_hess_l"];
            //     time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
            //     std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
            //     std::cout << "Total seconds from Ipopt w/o function: " << time_just_solver << std::endl;
            // }
            // else if (nonlinear_solver_name_ == "snopt")
            // {
            //     time_from_funcs += (double)stats["t_wall_nlp_grad"] + (double)stats["t_wall_nlp_jac_f"] + (double)stats["t_wall_nlp_jac_g"];
            //     time_just_solver += (double)stats["t_wall_total"] - time_from_funcs;
            //     std::cout << "Total seconds from Casadi functions: " << time_from_funcs << std::endl;
            //     std::cout << "Total seconds from SNOPT w/o function: " << time_just_solver << std::endl;
            // }
            prev_solution_ = result["x"];
            casadi::DM p_sol = arg["p"];

            for (size_t i = 0; i < trajectory_.size(); ++i)
            {
                casadi::DM seg_sol = prev_solution_(casadi::Slice(casadi_int(std::get<0>(ranges_decision_variables_[i])), casadi_int(std::get<1>(ranges_decision_variables_[i]))));
                casadi::DM seg_p_sol = p_sol(casadi::Slice(casadi_int(std::get<0>(ranges_parameters_[i])), casadi_int(std::get<1>(ranges_parameters_[i]))));
                casadi::DMVector sol_i = trajectory_[i]->ExtractSolution(seg_sol, seg_p_sol);
                if (i == 0)
                    sol_ = sol_i;
                else
                {
                    sol_[0] = horzcat(sol_[0], sol_i[0]);
                    sol_[1] = horzcat(sol_[1], sol_i[1]);
                }
            }
            return sol_;
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

                std::vector<double> state_times_vec = pseg->getSegmentTimes().get_elements();
                std::vector<double> input_times_vec = pseg->getUSegmentTimes().get_elements();
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