#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
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