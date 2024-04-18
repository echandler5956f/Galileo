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
            bool Solution::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result, GetSolutionType solution_type, AccessSolutionError &sol_error) const
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
                if (solution_type != STATE_ONLY && solution_type != INPUT_ONLY && solution_type != STATE_AND_INPUT)
                {
                    sol_error = AccessSolutionError::INVALID_SOLUTION_TYPE;
                    return false;
                }

                for (size_t i = 0; i < query_times.size(); i++)
                {
                    for (size_t j = 0; j < solution_segments_.size(); j++)
                    {
                        if ((query_times(i) > solution_segments_[j].initial_time) && (query_times(i) < solution_segments_[j].end_time))
                        {
                            if (solution_type == STATE_ONLY || solution_type == STATE_AND_INPUT)
                            {
                                Eigen::VectorXd times = solution_segments_[j].state_times;
                                int degree = solution_segments_[j].state_degree;
                                int deg = degree + 1;
                                size_t index = ((query_times(i) > times.array()).count() - 1) / deg;
                                Eigen::MatrixXd terms = solution_segments_[j].solx_segment.block(0, index * deg, solution_segments_[j].solx_segment.rows(), deg);
                                double knot_start_time = times[index * deg];
                                double knot_end_time = times[(index * deg) + deg - 1];
                                double scaled_time = (query_times(i) - knot_start_time) / (knot_end_time - knot_start_time);
                                state_result.col(i) = solution_segments_[j].state_poly.BarycentricInterpolation(scaled_time, terms);
                            }
                            if (solution_type == INPUT_ONLY || solution_type == STATE_AND_INPUT)
                            {
                                Eigen::VectorXd times = solution_segments_[j].input_times;
                                int degree = solution_segments_[j].input_degree;
                                int deg = degree + 1;
                                size_t index = ((query_times(i) > times.array()).count() - 1) / deg;
                                Eigen::MatrixXd terms = solution_segments_[j].solu_segment.block(0, index * deg, solution_segments_[j].solu_segment.rows(), deg);
                                double knot_start_time = times[index * deg];
                                double knot_end_time = times[(index * deg) + deg - 1];
                                double scaled_time = (query_times(i) - knot_start_time) / (knot_end_time - knot_start_time);
                                input_result.col(i) = solution_segments_[j].input_poly.BarycentricInterpolation(scaled_time, terms);
                            }
                            break;
                        }
                    }
                }
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
                    std::cout << "No segment decision times provided" << std::endl;
                    return false;
                }
                if (solution_segments_.size() == 0)
                {
                    sol_error = AccessSolutionError::SOLUTION_DNE;
                    std::cout << "No solution segments provided" << std::endl;
                    return false;
                }

                // Calculate the required size for w
                size_t w_size = 0;
                for (size_t i = 0; i < size_t(segment_decision_times.size() / 2); ++i)
                {
                    w_size += solution_segments_[i].solx_segment.rows() * segment_decision_times[2 * i].size1();
                    w_size += solution_segments_[i].solu_segment.rows() * segment_decision_times[2 * i + 1].size1();
                }

                // Resize w
                w.resize(w_size);
                Eigen::MatrixXd dummy_matrix;

                size_t index = 0;
                // #pragma omp parallel for
                for (size_t i = 0; i < size_t(segment_decision_times.size() / 2); ++i)
                {
                    casadi::DM x_phase_i_times = segment_decision_times[2 * i] + dt;
                    casadi::DM u_phase_i_times = segment_decision_times[2 * i + 1] + dt;

                    Eigen::VectorXd x_phase_i_times_eigen;
                    Eigen::VectorXd u_phase_i_times_eigen;

                    tools::casadiToEigen(x_phase_i_times, x_phase_i_times_eigen);
                    tools::casadiToEigen(u_phase_i_times, u_phase_i_times_eigen);

                    Eigen::MatrixXd segx_result(solution_segments_[i].solx_segment.rows(), x_phase_i_times_eigen.size());
                    GetSolution(x_phase_i_times_eigen, segx_result, dummy_matrix, STATE_ONLY, sol_error);

                    Eigen::MatrixXd segu_result(solution_segments_[i].solu_segment.rows(), u_phase_i_times_eigen.size());
                    GetSolution(u_phase_i_times_eigen, dummy_matrix, segu_result, INPUT_ONLY, sol_error);

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