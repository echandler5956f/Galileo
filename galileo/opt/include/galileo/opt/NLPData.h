#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
        struct NLPData
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
             * @brief Vector of general constraint lower bounds.
             *
             */
            std::vector<double> lbg;

            /**
             * @brief Vector of general constraint upper bounds.
             *
             */
            std::vector<double> ubg;

            /**
             * @brief Lower bounds associated with the decision variables.
             *
             */
            std::vector<double> lbw;

            /**
             * @brief Upper bounds associated with the decision variables.
             *
             */
            std::vector<double> ubw;

            /**
             * @brief Initial guess for the decision variables.
             */
            std::vector<double> w0;

            /**
             * @brief Expression for objective cost.
             *
             */
            casadi::MX J;
        };
    }
}