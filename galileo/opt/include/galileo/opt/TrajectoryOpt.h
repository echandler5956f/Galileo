#pragma once

#include "galileo/opt/PseudospectralSegment.h"
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"

namespace galileo
{
    namespace opt
    {
        /**
         * @brief The trajectory optimization class. This class is responsible for
            initializing the finite elements, and optimizing the trajectory.
         *
         */
        class TrajectoryOpt
        {
        public:
            /**
             * @brief Construct a new Trajectory Opt object.
             *
             * @param opts_ Options to pass to the solver
             * @param state_indices_ Helper class to get the state indices
             * @param problem Problem data containing the objective function and dynamics
             */
            TrajectoryOpt(casadi::Dict opts_, std::shared_ptr<States> state_indices_, std::shared_ptr<GeneralProblemData> problem);

            /**
             * @brief Initialize the finite elements.
             *
             * @param d The degree of the finite element polynomials
             * @param X0 The initial state to deviate from
             */
            void init_finite_elements(int d, casadi::DM X0);

            /**
             * @brief Optimize and return the solution.
             *
             * @return casadi::SXVector The solution
             */
            casadi::MXVector optimize();

            /**
             * @brief Get the times where the decision variables are evaluated.
            */
            std::vector<double> get_all_times() const;

        private:
            /**
             * @brief A Trajectory is made up of pseudospectral finite elements.
             *
             */
            std::vector<std::shared_ptr<PseudospectralSegment>> trajectory;

            /**
             * @brief Continuous-time function. The decision variables are infinitesimal deviations from the initial state,
                allowing for states to lie on a manifold. Fint is the function which maps these
                deviations back to the actual state space.
             *
             */
            casadi::Function Fint;

            /**
             * @brief Continuous-time function. The ineverse function of Fint. This is used to generate the initial guess for the states.
             *
             */
            casadi::Function Fdif;

            /**
             * @brief Continuous-time function. This function stores the system dynamics.
             *
             */
            casadi::Function F;

            /**
             * @brief The "running" or integrated cost function.
             *
             */
            casadi::Function L;

            /**
             * @brief The terminal cost function.
             *
             */
            casadi::Function Phi;

            /**
             * @brief Casadi solver options.
             *
             */
            casadi::Dict opts;

            /**
             * @brief Nonlinear function solver.
             *
             */
            casadi::Function solver;

            /**
             * @brief Slicer to get the states.
             *
             */
            std::shared_ptr<States> state_indices;

            /**
             * @brief Fixed time horizon of the entire trajectory.
             *
             */
            double T;

            /**
             * @brief Vector of all decision variables.
             *
             */
            casadi::MXVector w;

            /**
             * @brief Vector of all constraint expressions.
             *
             */
            casadi::MXVector g;

            /**
             * @brief Vector of all general constraint lower bounds.
             *
             */
            std::vector<double> lbg;

            /**
             * @brief Vector of all general constraint upper bounds.
             *
             */
            std::vector<double> ubg;

            /**
             * @brief  Lower bounds associated with the decision variables.
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

            /**
             * @brief Vector of all times where decision variables are evaluated.
             *
             */
            std::vector<double> all_times;

            /**
             * @brief Vector of plottable states
             *
             */
            casadi::MXVector xplot;

            /**
             * @brief Vector of plottable inputs
             *
             */
            casadi::MXVector uplot;
        };
    }
}