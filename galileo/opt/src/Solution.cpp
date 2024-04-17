#include "galileo/opt/Solution.h"

namespace galileo
{
    namespace opt
    {
        namespace solution
        {
            void Solution::UpdateSolution(std::vector<solution_segment_data_t> solution_segments)
            {
                solution_segments_ = solution_segments;
            }

            Eigen::MatrixXd Solution::GetSolutionSegment(const Eigen::VectorXd &query_times, const Eigen::MatrixXd &segment, const Eigen::VectorXd &times, const int degree, const LagrangePolynomial &poly) const
            {
                Eigen::MatrixXd result(segment.rows(), query_times.size());

                for (size_t i = 0; i < query_times.size(); i++)
                {
                    for (size_t j = 0; j < solution_segments_.size(); j++)
                    {
                        if (query_times(i) >= solution_segments_[j].initial_time && query_times(i) <= solution_segments_[j].end_time)
                        {
                            int deg = degree + 1;
                            size_t index = ((query_times(i) >= times.array()).count() - 1) / deg;
                            Eigen::MatrixXd terms = segment.block(0, index * deg, segment.rows(), deg);
                            double knot_start_time = times[index * deg];
                            double knot_end_time = times[(index * deg) + deg - 1];
                            double scaled_time = (query_times(i) - knot_start_time) / (knot_end_time - knot_start_time);
                            result.col(i) = poly.BarycentricInterpolation(scaled_time, terms);
                            break;
                        }
                    }
                }

                return result;
            }

            bool Solution::GetStateSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, AccessSolutionError &sol_error) const
            {
                if (query_times.size() == 0)
                {
                    sol_error = AccessSolutionError::NO_QUERY_TIMES_PROVIDED;
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    sol_error = AccessSolutionError::SOLUTION_DNE;
                    return false;
                }

                state_result = GetSolutionSegment(query_times, solution_segments_[0].solx_segment, solution_segments_[0].state_times, solution_segments_[0].state_degree, solution_segments_[0].state_poly);

                return true;
            }

            bool Solution::GetInputSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &input_result, AccessSolutionError &sol_error) const
            {
                if (query_times.size() == 0)
                {
                    sol_error = AccessSolutionError::NO_QUERY_TIMES_PROVIDED;
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    sol_error = AccessSolutionError::SOLUTION_DNE;
                    return false;
                }

                input_result = GetSolutionSegment(query_times, solution_segments_[0].solu_segment, solution_segments_[0].input_times, solution_segments_[0].input_degree, solution_segments_[0].input_poly);

                return true;
            }

            // state_result and input_result should be initialized to the correct size before calling GetSolution!
            bool Solution::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result, AccessSolutionError &sol_error) const
            {
                if (query_times.size() == 0)
                {
                    sol_error = AccessSolutionError::NO_QUERY_TIMES_PROVIDED;
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    sol_error = AccessSolutionError::SOLUTION_DNE;
                    return false;
                }

                state_result = GetSolutionSegment(query_times, solution_segments_[0].solx_segment, solution_segments_[0].state_times, solution_segments_[0].state_degree, solution_segments_[0].state_poly);
                input_result = GetSolutionSegment(query_times, solution_segments_[0].solu_segment, solution_segments_[0].input_times, solution_segments_[0].input_degree, solution_segments_[0].input_poly);

                return true;
            }

            /**
             * w is ordered as follows:
             *      [x_knots_phase0,
             *       x_collocation_points_phase0,
             *       u_knots_phase0,
             *       u_collocation_points_phase0,
             *       x_knots_phase1,
             *       ...,
             *       u_collocation_points_phaseN]
             *
             */
            bool Solution::GetNextGuessFromPrevSolution(const casadi::DMVector &segment_decision_times, Eigen::VectorXd &w, double dt, AccessSolutionError &sol_error) const
            {
                if (segment_decision_times.size() == 0)
                {
                    sol_error = AccessSolutionError::NO_QUERY_TIMES_PROVIDED;
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    sol_error = AccessSolutionError::SOLUTION_DNE;
                    return false;
                }

                // Calculate the required size for w
                size_t w_size = 0;
                for (size_t i = 0; i < segment_decision_times.size(); ++i)
                {
                    w_size += solution_segments_[i].solx_segment.rows() * segment_decision_times[2 * i].size1();
                    w_size += solution_segments_[i].solu_segment.rows() * segment_decision_times[2 * i + 1].size1();
                }

                // Resize w
                w.resize(w_size);

                size_t index = 0;
                #pragma omp parallel for
                for (size_t i = 0; i < segment_decision_times.size(); ++i)
                {
                    casadi::DM x_phase_i_times = segment_decision_times[2 * i] + dt;
                    casadi::DM u_phase_i_times = segment_decision_times[2 * i + 1] + dt;

                    Eigen::VectorXd x_phase_i_times_eigen;
                    Eigen::VectorXd u_phase_i_times_eigen;

                    tools::casadiToEigen(x_phase_i_times, x_phase_i_times_eigen);
                    tools::casadiToEigen(u_phase_i_times, u_phase_i_times_eigen);

                    Eigen::MatrixXd segx_result(solution_segments_[i].solx_segment.rows(), x_phase_i_times_eigen.size());
                    GetStateSolution(x_phase_i_times_eigen, segx_result, sol_error);

                    Eigen::MatrixXd segu_result(solution_segments_[i].solu_segment.rows(), u_phase_i_times_eigen.size());
                    GetInputSolution(u_phase_i_times_eigen, segu_result, sol_error);
                    
                    segx_result.resize(segx_result.size(), 1);
                    w.block(index, 0, segx_result.size(), 1) = segx_result;
                    index += segx_result.size();
                    segu_result.resize(segu_result.size(), 1);
                    w.block(index, 0, segu_result.size(), 1) = segu_result;
                    index += segu_result.size();
                }
                return true;
            }

            void Solution::UpdateConstraints(std::vector<std::vector<galileo::opt::ConstraintData>> constarint_data_segments)
            {
                constraint_data_segments_ = constarint_data_segments;
            }

            std::vector<std::vector<constraint_evaluations_t>> Solution::GetConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
            {
                GetSolution(query_times, state_result, input_result);

                std::vector<std::vector<constraint_evaluations_t>> constraint_evaluations;
                std::vector<constraint_evaluations_t> phase_constraint_evaluations;

                casadi::DM dm_state_result;
                casadi::DM dm_input_result;
                casadi::DM dm_times;

                tools::eigenToCasadi(state_result, dm_state_result);
                tools::eigenToCasadi(input_result, dm_input_result);
                tools::eigenToCasadi(query_times, dm_times);

                for (size_t i = 0; i < constraint_data_segments_.size(); ++i)
                {
                    phase_constraint_evaluations.clear();
                    std::vector<galileo::opt::ConstraintData> G = constraint_data_segments_[i];
                    tuple_size_t seg_range = getSegmentIndices(query_times, solution_segments_[i].initial_time, solution_segments_[i].end_time);
                    for (size_t j = 0; j < G.size(); ++j)
                    {
                        ConstraintData con_data = G[j];
                        casadi_int start_idx = casadi_int(std::get<0>(seg_range));
                        casadi_int end_idx = casadi_int(std::get<1>(seg_range));

                        casadi::DM con_eval = con_data.G.map(end_idx - start_idx)(casadi::DMVector{
                                                                                      dm_state_result(casadi::Slice(0, dm_state_result.rows()), casadi::Slice(start_idx, end_idx)),
                                                                                      dm_input_result(casadi::Slice(0, dm_input_result.rows()), casadi::Slice(start_idx, end_idx))})
                                                  .at(0);
                        casadi::DM con_lb = con_data.lower_bound.map(end_idx - start_idx)(casadi::DMVector{
                                                                                              dm_times(casadi::Slice(start_idx, end_idx), casadi::Slice(0, dm_times.columns()))})
                                                .at(0);
                        casadi::DM con_ub = con_data.upper_bound.map(end_idx - start_idx)(casadi::DMVector{
                                                                                              dm_times(casadi::Slice(start_idx, end_idx), casadi::Slice(0, dm_times.columns()))})
                                                .at(0);

                        Eigen::MatrixXd eval;
                        tools::casadiToEigen(con_eval, eval);
                        Eigen::MatrixXd lb;
                        tools::casadiToEigen(con_lb, lb);
                        Eigen::MatrixXd ub;
                        tools::casadiToEigen(con_ub, ub);

                        eval.transposeInPlace();
                        lb.transposeInPlace();
                        ub.transposeInPlace();

                        constraint_evaluations_t con_evals;
                        con_evals.metadata = con_data.metadata;
                        con_evals.times = query_times.block(std::get<0>(seg_range), 0, std::get<1>(seg_range) - std::get<0>(seg_range), 1);
                        con_evals.evaluation = eval;
                        con_evals.lower_bounds = lb;
                        con_evals.upper_bounds = ub;
                        phase_constraint_evaluations.push_back(con_evals);
                    }

                    constraint_evaluations.push_back(phase_constraint_evaluations);
                }

                return constraint_evaluations;
            }

            tuple_size_t Solution::getSegmentIndices(const Eigen::VectorXd &times, double start_time, double end_time) const
            {
                // size_t start_idx = (times.array() >= start_time).count() - 1;
                // size_t end_idx = (times.array() >= end_time).count() - 1;
                // return std::make_tuple(start_idx, end_idx);
                auto start_it = std::lower_bound(times.data(), times.data() + times.size(), start_time);
                auto end_it = std::upper_bound(times.data(), times.data() + times.size(), end_time);
                return std::make_tuple(std::distance(times.data(), start_it), std::distance(times.data(), end_it));
            }
        }
    }
}