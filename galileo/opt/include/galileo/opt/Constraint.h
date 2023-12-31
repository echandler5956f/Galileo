#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {

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
             * @param Fdif_ Continuous-time function. The ineverse function of Fint. This is used to generate the initial guess for the states.
             * @param F_ Dynamics of the system
             * @param L_ Running cost
             * @param Phi_ Terminal cost
             */
            GeneralProblemData(casadi::Function Fint_, casadi::Function Fdif_, casadi::Function F_, casadi::Function L_, casadi::Function Phi_)
            {
                this->Fint = Fint_;
                this->Fdif = Fdif_;
                this->F = F_;
                this->L = L_;
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
        };

        /**
         * @brief Results that describe a "built" constraint. 
         * This contains the constraint Function, Bounds, etcetera. 
         * 
         */
        struct ConstraintData
        {

            /**
             * @brief If global is true, use a map to apply the constraints across all knot points AND collocation points, and do not check flags.
             * 
             */
            bool global;

            /**
             * @brief If global is false, apply the constraints at these knot points (NOT collocation points) indices.
             * 
             */
            Eigen::VectorXi apply_at;

            /**
             * @brief Upper bounds of the constraint function.
             * 
             */
            casadi::Function upper_bound;

            /**
             * @brief Lower bounds of the constraint function.
             * 
             */
            casadi::Function lower_bound;

            /**
             * @brief The constraint function.
             * 
             */
            casadi::Function G;
        };

        /**
         * @brief Data for the decision variables.
         * 
         */
        struct DecisionData
        {

            /**
             * @brief Upper bound on decision variables.
             * 
             */
            casadi::Function upper_bound;

            /**
             * @brief Lower bound on decision variables.
             * 
             */
            casadi::Function lower_bound;

            /**
             * @brief Initial guess function for state and input as a function of time. Note that an inverse law for Fint will be used to generate the initial guess for the states.
             * 
             */
            casadi::Function initial_guess;

            /**
             * @brief Decision variables at one knot or collocation point (x, u)
             * 
             */
            casadi::SX w;
        };


        /**
         * @brief Extend this class to implement constraints.
         * 
         */
        template <class ProblemData>
        class ConstraintBuilder
        {
            /**
             * @brief Construct a new Constraint Builder object
             * 
             */
            ConstraintBuilder() {}

            /**
             * @brief Destroy the Constraint Builder object
             * 
             */
            virtual ~ConstraintBuilder() = default;

            /**
             * @brief Build constraint data for a given problem data.
             * 
             * @param problem_data Problem specific data
             * @param constraint_data Constraint specific data
             */
            void BuildConstraint(const ProblemData &problem_data, ConstraintData &constraint_data)
            {
                CreateApplyAt(problem_data, constraint_data.apply_at);
                CreateBounds(problem_data, constraint_data.upper_bound, constraint_data.lower_bound);
                CreateFunction(problem_data, constraint_data.G);
            }

            /**
             * @brief Generate flags for each knot point
             * 
             * @param problem_data Problem specific data
             * @param apply_at What knot points (indices) to apply the constraint at
             */
            virtual void CreateApplyAt(const ProblemData &problem_data, Eigen::VectorXi &apply_at) const = 0;

            /**
             * @brief Generate bounds for a vector of points
             * 
             * @param problem_data Problem specific data 
             * @param upper_bound Upper bound function to return
             * @param lower_bound Lower bound function to return
             */
            virtual void CreateBounds(const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound) const = 0;

            /**
             * @brief Generate a function to evaluate each point
             * 
             * @param problem_data Problem specific data
             * @param G The constraint function to return
             */
            virtual void CreateFunction(const ProblemData &problem_data, casadi::Function &G) const = 0;
        };
    };

};