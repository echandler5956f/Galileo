#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
        using tuple_size_t = std::tuple<size_t, size_t>;
        /**
         * @brief Struct to hold data needed to describe a Nonlinear Programming (NLP) problem.
         *
         */
        struct NLPProblemData
        {
            /**
             * @brief Vector of decision variables.
             *
             */
            casadi::MXVector w;

            /**
             * @brief Vector of constraint expressions.
             *
             */
            casadi::MXVector g;

            /**
             * @brief Vector of constant parameters.
             *
             */
            casadi::MXVector p;

            /**
             * @brief Expression for objective cost.
             *
             */
            casadi::MX J;
        };

        struct NLP
        {
            /**
             * @brief Symbolic data needed to describe the NLP problem.
             *
             */
            NLPProblemData nlp_prob_data;

            /**
             * @brief Ranges of the decision variables within each phase
             *
             */
            std::vector<tuple_size_t> ranges_decision_variables;

            /**
             * @brief Ranges of the parameters within each phase
             *
             */
            std::vector<tuple_size_t> ranges_parameters;

            /**
             * @brief The NLP solver function.
             *
             */
            casadi::Function nlp_solver;
        };

        /**
         * @brief Struct to hold data used as input to the NLP solver.
         *
         */
        struct NLPInputData
        {
            /**
             * @brief Lower bounds associated with the decision variables.
             *
             */
            casadi::DMVector lbw;

            /**
             * @brief Upper bounds associated with the decision variables.
             *
             */
            casadi::DMVector ubw;

            /**
             * @brief Initial guess for the decision variables.
             */
            casadi::DMVector w0;

            /**
             * @brief Vector of general constraint lower bounds.
             *
             */
            casadi::DMVector lbg;

            /**
             * @brief Vector of general constraint upper bounds.
             *
             */
            casadi::DMVector ubg;

            /**
             * @brief Vector of initial guess for the parameters.
             *
             */
            casadi::DMVector p0;
        };
    }
}