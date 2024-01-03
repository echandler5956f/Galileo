#pragma once

#include <Eigen/Core>
#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace variables
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
             * @param Fint_
             * @param F_ Dynamics of the system
             * @param L_ Running cost
             * @param Phi_ Terminal cost
             */
            GeneralProblemData(casadi::Function Fint_, casadi::Function F_, casadi::Function L_, casadi::Function Phi_)
            {
                this->Fint = Fint_;
                this->F = F_;
                this->L = L_;
                this->Phi = Phi_;
            }

            /**
             * @brief
             *
             */
            casadi::Function Fint;

            /**
             * @brief
             *
             */
            casadi::Function F;

            /**
             * @brief
             *
             */
            casadi::Function L;

            /**
             * @brief
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
             * @brief If global is true, use a map to apply the constraints across all knot points AND collocation points, and do not check flags
             *
             */
            bool global;

            /**
             * @brief If global is false, apply the constraints at these knot points (NOT collocation points) indices
             *
             */
            Eigen::VectorXi apply_at;

            /**
             * @brief
             *
             */
            casadi::Function upper_bound;

            /**
             * @brief
             *
             */
            casadi::Function lower_bound;

            /**
             * @brief
             *
             */
            casadi::Function G;
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
             * @brief
             *
             * @param problem_data
             * @param constraint_data
             */
            virtual void BuildConstraint(const ProblemData &problem_data, int knot_index, ConstraintData &constraint_data)
            {
                CreateApplyAt(problem_data, knot_index, constraint_data.apply_at);
                CreateBounds(problem_data, knot_index, constraint_data.upper_bound, constraint_data.lower_bound);
                CreateFunction(problem_data, knot_index, constraint_data.G);
            }

        private:
            /**
             * @brief Generate flags for each knot point
             *
             * @param problem_data)
             * @param apply_at
             */
            virtual void CreateApplyAt(const ProblemData &problem_data, int knot_index, Eigen::VectorXi &apply_at) const = 0;

            /**
             * @brief Generate bounds for a vector of points
             *
             * @param problem_data
             * @param upper_bound
             * @param lower_bound
             */
            virtual void CreateBounds(const ProblemData &problem_data, int knot_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const = 0;

            /**
             * @brief Generate a function to evaluate each point
             *
             * @param problem_data
             * @param G
             */
            virtual void CreateFunction(const ProblemData &problem_data, int knot_index, casadi::Function &G) const = 0;
        };
    };

};