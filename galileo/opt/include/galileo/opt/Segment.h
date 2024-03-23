#pragma once

#include "galileo/opt/States.h"
#include "galileo/opt/Constraint.h"
#include <cassert>

namespace galileo
{
    namespace opt
    {

        using tuple_size_t = std::tuple<size_t, size_t>;

        /**
         * @brief Base class for a segment used in TrajectoryOpt.
         *
         */
        class Segment
        {
        public:
            /**
             * @brief Destroy the Segment object
             *
             */
            virtual ~Segment() = default;

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return casadi::MX Solution values
             */
            virtual casadi::MXVector extractSolution(casadi::MX &w) const = 0;

            /**
             * @brief Get the initial state.
             *
             * @return casadi::MX The initial state
             */
            virtual casadi::MX getInitialState() const = 0;

            /**
             * @brief Get the initial state deviant.
             *
             * @return casadi::MX The initial state deviant
             */
            virtual casadi::MX getInitialStateDeviant() const = 0;

            /**
             * @brief Get the final state deviant.
             *
             * @return casadi::MX The final state deviant
             */
            virtual casadi::MX getFinalStateDeviant() const = 0;

            /**
             * @brief Get the actual final state.
             *
             * @return casadi::MX The final state.
             */
            virtual casadi::MX getFinalState() const = 0;

            /**
             * @brief Get the segment times vector.
             *
             * @return casadi::DM The segment times vector
             */
            virtual casadi::DM getSegmentTimes() const = 0;

            /**
             * @brief Get the input segment times vector.
             * 
             * @return casadi::DM The input segment times vector
             */
            virtual casadi::DM getUSegmentTimes() const = 0;

            /**
             * @brief Initialize the vector of segment times which constraints are evaluated at.
             *
             * @param global_times Vector of global times
             */
            virtual void initializeSegmentTimeVector(casadi::DM &global_times) = 0;

            /**
             * @brief Initialize the vector of times which coincide to the decision variables U occur at.
             *
             */
            virtual void initializeInputTimeVector(casadi::DM &global_times) = 0;

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             * @param x0_local Local starting state to integrate from
             *
             */
            virtual void initializeKnotSegments(casadi::DM x0_global, casadi::MX x0_local) = 0;

            /**
             * @brief Build the function graph.
             *
             * @param G Vector of constraint data
             * @param Wdata Decision bound and initial guess data for the state and input
             */
            virtual void initializeExpressionGraph(std::vector<ConstraintData> G, std::shared_ptr<DecisionData> Wdata) = 0;

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param w Decision variable vector to fill
             * @param g Constraint vector to fill
             */
            virtual void evaluateExpressionGraph(casadi::MX &J0, casadi::MXVector &w, casadi::MXVector &g) = 0;

            /**
             * @brief Fills the lower bounds on decision variable (lbw) and upper bounds on decision variable (ubw) vectors with values.
             *
             * This function takes in two vectors, lbw and ubw, and fills them with values.
             * The filled values represent the constraint lower and upper bounds on decision variables.
             *
             * @param lbw The vector to be filled with lower bound values on decision variables.
             * @param ubw The vector to be filled with upper bound values on decision variables.
             */
            virtual void fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw) = 0;

            /**
             * @brief Fills the lower bounds on general constraints (lbg) and upper bounds on general constraints (ubg) vectors with values.
             *
             * This function takes in two vectors, lbg and ubg, and fills them with values.
             * The filled values represent the general constraint lower and upper bounds.
             *
             * @param lbg The vector to be filled with general lower bound values.
             * @param ubg The vector to be filled with general upper bound values.
             */
            virtual void fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg) = 0;

            /**
             * @brief Fills the initial guess vector (w0) with values.
             *
             * This function takes in a vector, w0, and fills it with values.
             * The filled values represent the initial guess for the decision variables.
             *
             * @param w0 The vector to be filled with initial guess values.
             */
            virtual void fill_w0(std::vector<double> &w0) const = 0;

            /**
             * @brief Returns the starting and ending index in w.
             *
             * @return tuple_size_t The range of indices
             */
            virtual tuple_size_t get_range_idx_decision_variables() const = 0;

            /**
             * @brief Returns the starting and ending index in g (call after evaluate_expression_graph!).
             *
             * @return tuple_size_t The range of indices
             */
            virtual tuple_size_t get_range_idx_constraint_expressions() const = 0;

            /**
             * @brief Returns the starting and ending index in lbg/ubg.
             *
             * @return tuple_size_t The range of indices
             */
            virtual tuple_size_t get_range_idx_constraint_bounds() const = 0;

            /**
             * @brief Returns the starting and ending index in lbw/ubw.
             *
             * @return tuple_size_t The range of indices
             */
            virtual tuple_size_t get_range_idx_decision_bounds() const = 0;

        public:
            /**
             * @brief Local initial state.
             *
             */
            casadi::MX x0_local;

            /**
             * @brief Global initial state.
             *
             */
            casadi::DM x0_global;

            /**
             * @brief Period of EACH KNOT SEGMENT within this pseudospectral segment.
             *
             */
            double h;

            /**
             * @brief Total period (helper variable calculated from h and knot_num).
             *
             */
            double T;

            /**
             * @brief Starting and ending index of the decision variables in w corresponding to this segment.
             *
             */
            tuple_size_t w_range;

            /**
             * @brief Starting and ending index of the constraint expressions in g corresponding to this segment.
             *
             */
            tuple_size_t g_range;

            /**
             * @brief Starting and ending index of the bounds in lbg/ubg corresponding to this segment.
             *
             */
            tuple_size_t lbg_ubg_range;

            /**
             * @brief Starting and ending index of the bounds in lbw/ubw corresponding to this segment.
             *
             */
            tuple_size_t lbw_ubw_range;
        };
    }
}