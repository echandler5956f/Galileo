#pragma once

#include <Eigen/Core>
#include <pinocchio/autodiff/casadi.hpp>

namespace acro
{
    namespace variables
    {

        /**
         * @brief 
         * 
         */
        struct ProblemData
        {
            /**
             * @brief Construct a new Problem Data object
             * 
             * @param Fint_ 
             * @param F_ 
             * @param L_ 
             * @param Phi_ 
             */
            ProblemData(casadi::Function Fint_, casadi::Function F_, casadi::Function L_, casadi::Function Phi_)
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
         * @brief 
         * 
         */
        struct ConstraintData
        {
            /**
             * @brief 
             * 
             */
            Eigen::VectorXi flags;

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
            casadi::Function F;
        };

        /**
         * @brief 
         * 
         */
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
            void BuildConstraint(const ProblemData &problem_data, ConstraintData &constraint_data)
            {
                CreateFlags(problem_data, constraint_data.flags);
                CreateBounds(problem_data, constraint_data.upper_bound, constraint_data.lower_bound);
                CreateFunction(problem_data, constraint_data.F);
            }

            /**
             * @brief Generate flags for each collocation point
             * 
             * @param problem_data 
             * @param flags 
             */
            virtual void CreateFlags(const ProblemData &problem_data, Eigen::VectorXi &flags) const;

            /**
             * @brief Generate bounds for a vector of concatinated collocation points
             * 
             * @param problem_data 
             * @param upper_bound 
             * @param lower_bound 
             */
            virtual void CreateBounds(const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound) const;

            /**
             * @brief Generate a function to evaluate each collocation point.
             * 
             * @param problem_data 
             * @param F 
             */
            virtual void CreateFunction(const ProblemData &problem_data, casadi::Function &F) const;
        };
    };

};