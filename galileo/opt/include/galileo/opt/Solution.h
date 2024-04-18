#pragma once

#include "galileo/opt/LagrangePolynomial.h"
#include "galileo/opt/ProblemData.h"
#include "galileo/tools/CasadiTools.h"
#include "galileo/opt/PhaseSequence.h"
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

                enum GetSolutionType
                {
                    STATE_ONLY,
                    INPUT_ONLY,
                    STATE_AND_INPUT
                };

                enum AccessSolutionError
                {
                    OK,
                    NO_QUERY_TIMES_PROVIDED,
                    SOLUTION_DNE,
                    INVALID_SOLUTION_TYPE
                };

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
                    GetSolutionType solution_type = STATE_AND_INPUT;
                    return GetSolution(query_times, state_result, input_result, solution_type, sol_error);
                }

                /**
                 * @brief Get the solution at a set of query times.
                 *
                 * @param query_times A vector of times at which to query the solution.
                 * @param state_result The state result at each query time.
                 * @param input_result The input result at each query time.
                 * @param solution_type The type of solution to get.
                 * @param sol_error An error code that is set if this fails to get a solution.
                 *
                 * @return bool True if the solution exists at the query times.
                 */
                bool GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result, GetSolutionType solution_type, AccessSolutionError &sol_error) const;

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
                 * @param query_times A vector of times at which to query the constraints.
                 * @param state_result The state at each query time.
                 * @param input_result The input at each query time.
                 *
                 * @return std::vector<std::vector<constraint_evaluations_t>> A vector of constraint evaluations.
                 */
                template <typename MODE_T>
                    std::vector<std::vector<constraint_evaluations_t>> GetConstraints(const Eigen::VectorXd &query_times, const Eigen::MatrixXd &state_result, const Eigen::MatrixXd &input_result, const std::shared_ptr<PhaseSequence<MODE_T>> sequence) const
                    {
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
                            tuple_size_t seg_range = getSegmentIndices(query_times, sequence->getPhaseStartTimes()[i], sequence->getPhaseStartTimes()[i] + sequence->getPhase(i).time_value);
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