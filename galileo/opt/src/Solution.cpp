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

            // state_result and input_result should be initialized to the correct size before calling GetSolution!
            bool Solution::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
            {
                if (query_times.size() == 0)
                {
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    return false;
                }

                auto clock_start_time = std::chrono::high_resolution_clock::now();

                for (size_t i = 0; i < query_times.size(); i++)
                {
                    for (size_t j = 0; j < solution_segments_.size(); j++)
                    {
                        if (query_times(i) >= solution_segments_[j].initial_time && query_times(i) <= solution_segments_[j].end_time)
                        {
                            int state_deg = solution_segments_[j].state_degree + 1;
                            size_t state_index = ((query_times(i) > solution_segments_[j].state_times.array()).count() - 1) / state_deg;
                            Eigen::MatrixXd state_terms = solution_segments_[j].solx_segment.block(0, state_index * state_deg, solution_segments_[j].solx_segment.rows(), state_deg);
                            double state_knot_start_time = solution_segments_[j].state_times[state_index * state_deg - 1];
                            double state_knot_end_time = solution_segments_[j].state_times[(state_index * state_deg) + state_deg - 1];

                            double state_scaled_time = (query_times(i) - state_knot_start_time) / (state_knot_end_time - state_knot_start_time);
                            state_result.col(i) = solution_segments_[j].state_poly.barycentricInterpolation(state_scaled_time, state_terms);

                            int input_deg = solution_segments_[j].input_degree + 1;
                            size_t input_index = ((query_times(i) > solution_segments_[j].input_times.array()).count() - 1) / input_deg;
                            Eigen::MatrixXd input_terms = solution_segments_[j].solu_segment.block(0, input_index * state_deg, solution_segments_[j].solu_segment.rows(), input_deg);
                            double input_knot_start_time = solution_segments_[j].input_times[input_index * input_deg - 1];
                            double input_knot_end_time = solution_segments_[j].input_times[(input_index * input_deg) + input_deg - 1];
                            double input_scaled_time = (query_times(i) - input_knot_start_time) / (input_knot_end_time - input_knot_start_time);
                            input_result.col(i) = solution_segments_[j].input_poly.barycentricInterpolation(input_scaled_time, input_terms);
                            break;
                        }
                    }
                }

                auto clock_end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_time = clock_end_time - clock_start_time;
                std::cout << "GetSolution took " << elapsed_time.count() << " seconds" << std::endl;
                return true;
            }

            void Solution::UpdateConstraints(std::vector<std::vector<galileo::opt::ConstraintData>> constarint_data_segments)
            {
                constraint_data_segments_ = constarint_data_segments;
            }

            // specialized implementation using memcpy; not really needed elsewhere
            void eigenToCasadi(const Eigen::MatrixXd &eigen_matrix, casadi::DM &casadi_matrix)
            {
                size_t rows = eigen_matrix.rows();
                size_t cols = eigen_matrix.cols();

                casadi_matrix.resize(rows, cols);
                casadi_matrix = casadi::DM::zeros(rows, cols);

                std::memcpy(casadi_matrix.ptr(), eigen_matrix.data(), sizeof(double) * rows * cols);
            }

            // specialized implementation using memcpy; not really needed elsewhere
            void casadiToEigen(const casadi::DM &casadi_matrix, Eigen::MatrixXd &eigen_matrix)
            {
                size_t rows = casadi_matrix.size1();
                size_t cols = casadi_matrix.size2();

                eigen_matrix.resize(rows, cols);
                eigen_matrix.setZero(rows, cols);

                std::memcpy(eigen_matrix.data(), casadi_matrix.ptr(), sizeof(double) * rows * cols);
            }

            std::vector<std::vector<constraint_evaluations_t>> Solution::GetConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
            {
                GetSolution(query_times, state_result, input_result);
                auto clock_start_time = std::chrono::high_resolution_clock::now();
                std::vector<std::vector<constraint_evaluations_t>> constraint_evaluations;
                std::vector<constraint_evaluations_t> phase_constraint_evaluations;

                casadi::DM dm_state_result;
                casadi::DM dm_input_result;
                casadi::DM dm_times;

                eigenToCasadi(state_result, dm_state_result);
                eigenToCasadi(input_result, dm_input_result);
                eigenToCasadi(query_times, dm_times);

                std::cout << "dm_state_result: " << dm_state_result.size1() << " x " << dm_state_result.size2() << std::endl;
                std::cout << "dm_input_result: " << dm_input_result.size1() << " x " << dm_input_result.size2() << std::endl;
                std::cout << "dm_times: " << dm_times.size1() << " x " << dm_times.size2() << std::endl;

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
                        casadiToEigen(con_eval, eval);
                        Eigen::MatrixXd lb;
                        casadiToEigen(con_lb, lb);
                        Eigen::MatrixXd ub;
                        casadiToEigen(con_ub, ub);

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
                auto clock_end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_time = clock_end_time - clock_start_time;
                std::cout << "GetConstraints took " << elapsed_time.count() << " seconds" << std::endl;
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