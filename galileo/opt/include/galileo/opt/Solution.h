#pragma once

#include "galileo/opt/LagrangePolynomial.h"
#include <Eigen/Dense>
#include <string>

#include <chrono>

namespace galileo
{
    namespace opt
    {
        namespace solution
        {
            using tuple_size_t = std::tuple<size_t, size_t>;

            struct constraint_metadata_t
            {
                /**
                 * @brief Name of the constraint.
                 *
                 */
                std::string name;

                /**
                 * @brief A vector of index ranges which group certain constraint rows together to make plotting easier. An empty vector means to include all constraint rows in the same plot.
                 *
                 */
                std::vector<tuple_size_t> plot_groupings;

                /**
                 * @brief Titles of the plots.
                 *
                 */
                std::vector<std::string> plot_titles;

                /**
                 * @brief Names of the plots on the legend.
                 *
                 */
                std::vector<std::vector<std::string>> plot_names;
            };

            struct constraint_evaluations_t
            {
                Eigen::VectorXd times;
                Eigen::MatrixXd evaluation;
                Eigen::MatrixXd lower_bounds;
                Eigen::MatrixXd upper_bounds;
                constraint_metadata_t metadata;
            };

            struct solution_t
            {
                solution_t() {}
                solution_t(Eigen::VectorXd times_) { this->times = times_; }
                solution_t(Eigen::VectorXd times_, Eigen::MatrixXd state_result_) { this->times = times_; this->state_result = state_result_; }
                solution_t(Eigen::VectorXd times_, Eigen::MatrixXd state_result_, Eigen::MatrixXd input_result_) { this->times = times_; this->state_result = state_result_; this->input_result = input_result_; }
                Eigen::VectorXd times;
                Eigen::MatrixXd state_result;
                Eigen::MatrixXd input_result;
            };

            struct solution_segment_data_t
            {
                Eigen::VectorXd state_times;
                Eigen::MatrixXd solx_segment;
                int state_degree;

                Eigen::VectorXd input_times;
                Eigen::MatrixXd solu_segment;
                int input_degree;

                double initial_time;
                double end_time;
                int num_knots;

                LagrangePolynomial state_poly;
                LagrangePolynomial input_poly;
            };

            class Solution
            {
            public:
                Solution(){};

                void UpdateSolution(std::vector<solution_segment_data_t> solution_segments);

                void GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const;

            // protected:
                // Eigen::VectorXd current_state_times_;
                // Eigen::MatrixXd current_state_sol_;

                // Eigen::VectorXd current_input_times_;
                // Eigen::MatrixXd current_input_sol_;

                std::vector<solution_segment_data_t> solution_segments_;
            };
        }
    }
}