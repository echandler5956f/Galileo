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
            void Solution::GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
            {
                std::cout << "GetSolution called" << std::endl;
                auto clock_start_time = std::chrono::high_resolution_clock::now();

                for (int i = 0; i < query_times.size(); i++)
                {
                    for (int j = 0; j < solution_segments_.size(); j++)
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
            }
        }
    }
}