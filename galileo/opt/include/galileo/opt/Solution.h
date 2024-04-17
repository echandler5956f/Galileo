#pragma once

#include "galileo/opt/LagrangePolynomial.h"
#include "galileo/opt/ProblemData.h"
#include "galileo/tools/CasadiTools.h"
#include <Eigen/Dense>
#include <string>

#include <omp.h>

#include <chrono>

namespace galileo
{
    namespace opt
    {
        namespace solution
        {
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
                 * @param query_times The times at which to query the solution.
                 * @param segment The segment to query.
                 * @param times The times at which the solution is evaluated.
                 * @param degree The degree of the polynomial.
                 * @param poly The polynomial to use.
                 * @return Eigen::MatrixXd The solution at the query times.
                 */
                Eigen::MatrixXd GetSolutionSegment(const Eigen::VectorXd &query_times, const Eigen::MatrixXd &segment, const Eigen::VectorXd &times, const int degree, const LagrangePolynomial &poly) const;

                enum AccessSolutionError
                {
                    OK,
                    NO_QUERY_TIMES_PROVIDED,
                    SOLUTION_DNE
                };

                /**
                 * @brief Get the state solution at a set of query times.
                 *
                 * @param query_times The times at which to query the solution.
                 * @param state_result The state result at each query time.
                 * @param sol_error An error code that is set if this fails to get a solution.
                 * @return true Success
                 * @return false Failure
                 */
                bool GetStateSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, AccessSolutionError &sol_error) const;

                /**
                 * @brief Get the input solution at a set of query times.
                 *
                 * @param query_times The times at which to query the solution.
                 * @param input_result The input result at each query time.
                 * @param sol_error An error code that is set if this fails to get a solution.
                 * @return true Success
                 * @return false Failure
                 */
                bool GetInputSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &input_result, AccessSolutionError &sol_error) const;

                /**
                 * @brief Get the solution at a set of query times.
                 *
                 * @param query_times A vector of times at which to query the solution.
                 * @param state_result The state result at each query time.
                 * @param input_result The input result at each query time.
                 *
                 * @return bool True if the solution exists at the query times.
                 */
                bool GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const
                {
                    AccessSolutionError sol_error;
                    return GetSolution(query_times, state_result, input_result, sol_error);
                }

                /**
                 * @brief Get the solution at a set of query times.
                 *
                 * @param query_times A vector of times at which to query the solution.
                 * @param state_result The state result at each query time.
                 * @param input_result The input result at each query time.
                 * @param sol_error An error code that is set if this fails to get a solution.
                 *
                 * @return bool True if the solution exists at the query times.
                 */
                bool GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result, AccessSolutionError &sol_error) const;

                /**
                 * @brief Get the next guess from the previous solution.
                 *
                 * @param segment_decision_times The decision times of the segments.
                 * @param w The next guess.
                 * @param dt The time step.
                 * @return true Success
                 * @return false Failure
                 */
                bool GetNextGuessFromPrevSolution(const casadi::DMVector &segment_decision_times, Eigen::VectorXd &w, double dt) const
                {
                    AccessSolutionError sol_error;
                    return GetNextGuessFromPrevSolution(segment_decision_times, w, dt, sol_error);
                }

                /**
                 * @brief Get the next guess from the previous solution.
                 *
                 * @param segment_decision_times The decision times of the segments.
                 * @param w The next guess.
                 * @param dt The time step.
                 * @param sol_error An error code that is set if this fails to get a solution.
                 * @return true Success
                 * @return false Failure
                 */
                bool GetNextGuessFromPrevSolution(const casadi::DMVector &segment_decision_times, Eigen::VectorXd &w, double dt, AccessSolutionError &sol_error) const;

                /**
                 * @brief Update the constraints with new constraint data segments.
                 *
                 * @param constarint_data_segments A vector of constraint data segments.
                 */
                void UpdateConstraints(std::vector<std::vector<galileo::opt::ConstraintData>> constarint_data_segments);

                /**
                 * @brief Get the constraint evaluations at a set of query times.
                 *
                 * You should call either GetSolution or GetConstraints, but not both.
                 * Calling GetConstraints will call GetSolution internally and fill in the state and input results.
                 *
                 * @param query_times A vector of times at which to query the constraints.
                 * @param state_result The state result at each query time. This is found by calling GetSolution.
                 * @param input_result The input result at each query time. This is found by calling GetSolution.
                 *
                 * @return std::vector<std::vector<constraint_evaluations_t>> A vector of constraint evaluations.
                 */
                std::vector<std::vector<constraint_evaluations_t>> GetConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result) const;

                /**
                 * @brief Get the segment indices for a given time range.
                 *
                 * @param times The times at which the solution is evaluated.
                 * @param start_time Start time of the segment
                 * @param end_time End time of the segment
                 * @return tuple_size_t A tuple of the start and end indices of the segment in the times vector.
                 */
                tuple_size_t getSegmentIndices(const Eigen::VectorXd &times, double start_time, double end_time) const;

                /**
                 * @brief Get if this contains a solution
                 *
                 * @return bool True if the solution set is not empty.
                 */
                bool isSolutionSet() const { return !solution_segments_.empty(); }

            private:
                /**
                 * @brief The solution segments.
                 *
                 */
                std::vector<solution_segment_data_t> solution_segments_;

                /**
                 * @brief The constraint data segments.
                 *
                 */
                std::vector<std::vector<galileo::opt::ConstraintData>> constraint_data_segments_;
            };
        }
    }
}