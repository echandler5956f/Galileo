#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Struct to hold data for a Nonlinear Programming (NLP) problem.
         * 
         */
        struct NLPData
        {
            /**
             * @brief Vector of decision variables.
             *
             */
            casadi::MXVector w;

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
             * @brief Vector of constraint expressions.
             *
             */
            casadi::MXVector g;

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
             * @brief Expression for objective cost.
             *
             */
            casadi::MX J;
        };
    }
}