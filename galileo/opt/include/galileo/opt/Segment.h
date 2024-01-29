#pragma once

#include "galileo/opt/States.h"
#include "galileo/opt/Constraint.h"
#include <vector>
#include <string>
#include <cassert>
#include <memory>

using namespace casadi;

namespace galileo
{
    namespace opt
    {

        using tuple_size_t = std::tuple<std::size_t, std::size_t>;

        class Segment
        {
        public:
            virtual ~Segment() = default;

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return MX Solution values
             */
            virtual MXVector extract_solution(MX &w) const = 0;

            /**
             * @brief Get the initial state.
             *
             * @return MX The initial state
             */
            virtual MX get_initial_state() const = 0;

            /**
             * @brief Get the initial state deviant.
             *
             * @return MX The initial state deviant
             */
            virtual MX get_initial_state_deviant() const = 0;

            /**
             * @brief Get the final state deviant.
             *
             * @return MX The final state deviant
             */
            virtual MX get_final_state_deviant() const = 0;

            /**
             * @brief Get the actual final state.
             *
             * @return MX The final state.
             */
            virtual MX get_final_state() const = 0;

            /**
             * @brief Get the local times vector.
             *
             * @return casadi::DM The local times vector
             */
            virtual casadi::DM get_local_times() const = 0;

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
            MX x0_local;

            /**
             * @brief Global initial state.
             *
             */
            DM x0_global;

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