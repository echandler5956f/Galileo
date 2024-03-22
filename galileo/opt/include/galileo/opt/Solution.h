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

            /**
             * @brief Struct for storing metadata about a constraint.
             * 
             */
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

            /**
             * @brief Struct for storing constraint evaluations.
             * 
             */
            struct constraint_evaluations_t
            {
                /**
                 * @brief A vector of times at which the constraint was evaluated.
                 *
                 */
                Eigen::VectorXd times;

                /**
                 * @brief The evaluation of the constraint at each time.
                 *
                 */
                Eigen::MatrixXd evaluation;

                /**
                 * @brief The lower bounds of the constraint at each time.
                 *
                 */
                Eigen::MatrixXd lower_bounds;

                /**
                 * @brief The upper bounds of the constraint at each time.
                 *
                 */
                Eigen::MatrixXd upper_bounds;

                /**
                 * @brief Plotting metadata for the constraint.
                 *
                 */
                constraint_metadata_t metadata;
            };

            /**
             * @brief Struct for storing a solution.
             * 
             */
            struct solution_t
            {
                /**
                 * @brief Default constructor.
                 *
                 */
                solution_t() {}

                /**
                 * @brief Construct a new solution_t object.
                 *
                 * @param times_ A vector of times at which the solution is evaluated.
                 */
                solution_t(Eigen::VectorXd times_) { this->times = times_; }

                /**
                 * @brief Construct a new solution_t object.
                 *
                 * @param times_ A vector of times at which the solution is evaluated.
                 * @param state_result_ The state result at each time.
                 */
                solution_t(Eigen::VectorXd times_, Eigen::MatrixXd state_result_)
                {
                    this->times = times_;
                    this->state_result = state_result_;
                }

                /**
                 * @brief Construct a new solution_t object.
                 *
                 * @param times_ A vector of times at which the solution is evaluated.
                 * @param state_result_ The state result at each time.
                 * @param input_result_ The input result at each time.
                 */
                solution_t(Eigen::VectorXd times_, Eigen::MatrixXd state_result_, Eigen::MatrixXd input_result_)
                {
                    this->times = times_;
                    this->state_result = state_result_;
                    this->input_result = input_result_;
                }

                /**
                 * @brief A vector of times at which the solution is evaluated.
                 *
                 */
                Eigen::VectorXd times;

                /**
                 * @brief The state result at each time.
                 *
                 */
                Eigen::MatrixXd state_result;

                /**
                 * @brief The input result at each time.
                 *
                 */
                Eigen::MatrixXd input_result;
            };

            /**
             * @brief Struct for storing solution segment data.
             * 
             */
            struct solution_segment_data_t
            {
                /**
                 * @brief A vector of times at which the state is evaluated for this segment.
                 *
                 */
                Eigen::VectorXd state_times;

                /**
                 * @brief The state solution for this segment.
                 *
                 */
                Eigen::MatrixXd solx_segment;

                /**
                 * @brief The degree of the state polynomial.
                 *
                 */
                int state_degree;

                /**
                 * @brief A vector of times at which the input is evaluated for this segment.
                 *
                 */
                Eigen::VectorXd input_times;

                /**
                 * @brief The input solution for this segment.
                 *
                 */
                Eigen::MatrixXd solu_segment;

                /**
                 * @brief The degree of the input polynomial.
                 *
                 */
                int input_degree;

                /**
                 * @brief The initial time of the segment.
                 *
                 */
                double initial_time;

                /**
                 * @brief The end time of the segment.
                 *
                 */
                double end_time;

                /**
                 * @brief The number of knots in the segment.
                 *
                 */
                int num_knots;

                /**
                 * @brief The state polynomial.
                 *
                 */
                LagrangePolynomial state_poly;

                /**
                 * @brief The input polynomial.
                 *
                 */
                LagrangePolynomial input_poly;
            };

            /**
             * @brief Class for storing and retrieving solutions.
             *
             */
            class Solution
            {
            public:
                /**
                 * @brief Default constructor.
                 *
                 */
                Solution(){};

                /**
                 * @brief Update the solution with new segments.
                 *
                 * @param solution_segments A vector of solution segments.
                 */
                void UpdateSolution(std::vector<solution_segment_data_t> solution_segments);

                /**
                 * @brief Get the solution at a set of query times.
                 *
                 * @param query_times A vector of times at which to query the solution.
                 * @param state_result The state result at each query time.
                 * @param input_result The input result at each query time.
                 */
                void GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const;

                /**
                 * @brief The solution segments.
                 *
                 */
                std::vector<solution_segment_data_t> solution_segments_;
            };
        }
    }
}