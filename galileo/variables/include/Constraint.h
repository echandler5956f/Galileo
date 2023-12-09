#pragma once

#include <Eigen/Core>
#include <pinocchio/autodiff/casadi.hpp>

namespace acro
{
    namespace variables
    {

        struct ProblemData
        {
            ProblemData(casadi::Function Fint_, casadi::Function F_, casadi::Function L_, casadi::Function Phi_)
            {
                this->Fint = Fint_;
                this->F = F_;
                this->L = L_;
                this->Phi = Phi_;
            }
            casadi::Function Fint;
            casadi::Function F;
            casadi::Function L;
            casadi::Function Phi;
        };

        struct ConstraintData
        {
            Eigen::VectorXi flags;

            casadi::Function upper_bound;
            casadi::Function lower_bound;

            casadi::Function F;
        };

        class ConstraintBuilder
        {
            ConstraintBuilder() {}

            virtual ~ConstraintBuilder() = default;

            void BuildConstraint(const ProblemData &problem_data, ConstraintData &constraint_data)
            {
                CreateFlags(problem_data, constraint_data.flags);
                CreateBounds(problem_data, constraint_data.upper_bound, constraint_data.lower_bound);
                CreateFunction(problem_data, constraint_data.F);
            }

            // Generate flags for each collocation point
            virtual void CreateFlags(const ProblemData &problem_data, Eigen::VectorXi &flags) const;

            // Generate bounds for a vector of concatinated collocation points
            virtual void CreateBounds(const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound) const;

            // Generate a function to evaluate each collocation point.
            virtual void CreateFunction(const ProblemData &problem_data, casadi::Function &F) const;
        };
    };

};