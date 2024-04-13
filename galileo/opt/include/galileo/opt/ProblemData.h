#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
        using tuple_size_t = std::tuple<size_t, size_t>;
        /**
         * @brief Data necessary to build the problem or constraints.
         * It is good practice to have a corresponding constraint specific
         * Problem Data for any data needed to build a constraint.
         *
         */
        struct GeneralProblemData
        {
            /**
             * @brief Construct a new Problem Data object.
             *
             * @param Fint_ Continuous-time function. The decision variables are infinitesimal deviations from the initial state,
                allowing for states to lie on a manifold. Fint is the function which maps these
                deviations back to the actual state space
             * @param Fdiff_ Continuous-time function. The ineverse function of Fint. This is used to generate the initial guess for the states.
             * @param Phi_ Terminal cost
             */
            GeneralProblemData(casadi::Function Fint_, casadi::Function Fdiff_, casadi::Function Phi_)
            {
                this->Fint = Fint_;
                this->Fdiff = Fdiff_;
                this->Phi = Phi_;
            }

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
            casadi::Function Fdiff;

            /**
             * @brief The terminal cost function.
             *
             */
            casadi::Function Phi;
        };

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
         * @brief Results that describe a "built" constraint.
         * This contains the constraint Function, Bounds, etc.
         *
         */
        struct ConstraintData
        {
            /**
             * @brief Lower bounds of the constraint function.
             *
             */
            casadi::Function lower_bound;

            /**
             * @brief Upper bounds of the constraint function.
             *
             */
            casadi::Function upper_bound;

            /**
             * @brief Constraint function.
             *
             */
            casadi::Function G;

            /**
             * @brief Metadata for the constraint.
             *
             */
            constraint_metadata_t metadata;
        };

        /**
         * @brief Data for the decision variables.
         *
         */
        struct DecisionData
        {
            /**
             * @brief Lower bound on decision variables.
             *
             */
            casadi::Function lower_bound;

            /**
             * @brief Upper bound on decision variables.
             *
             */
            casadi::Function upper_bound;

            /**
             * @brief Initial guess function for state and input as a function of time. Note that an inverse law for Fint will be used to generate the initial guess for the states.
             *
             */
            casadi::Function initial_guess;
        };

        /**
         * @brief Extend this class to implement constraints.
         *
         */
        template <class ProblemData>
        class ConstraintBuilder
        {
        public:
            /**
             * @brief Construct a new Constraint Builder object.
             *
             */
            ConstraintBuilder() {}

            /**
             * @brief Destroy the Constraint Builder object.
             *
             */
            virtual ~ConstraintBuilder() = default;

            /**
             * @brief Build constraint data for a given problem data.
             *
             * @param problem_data Problem specific data
             * @param phase_index Index to build constraint data for
             * @param constraint_data Constraint specific data
             */
            virtual void buildConstraint(const ProblemData &problem_data, int phase_index, ConstraintData &constraint_data) = 0;
        };

        /**
         * @brief Extend this class to implement constraints.
         *
         */
        template <class ProblemData>
        class DecisionDataBuilder
        {
        public:
            /**
             * @brief Construct a new Decision Data Builder object.
             *
             */
            DecisionDataBuilder() {}

            /**
             * @brief Destroy the Decision Data Builder object.
             *
             */
            virtual ~DecisionDataBuilder() = default;

            /**
             * @brief Build constraint data for a given problem data.
             *
             * @param problem_data Problem specific data
             * @param phase_index Index to build constraint data for
             * @param decision_data Constraint specific data
             */
            virtual void buildDecisionData(const ProblemData &problem_data, int phase_index, DecisionData &decision_data) = 0;
        };
    };
};