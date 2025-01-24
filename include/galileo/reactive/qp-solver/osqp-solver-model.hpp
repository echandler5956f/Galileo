#ifndef __galileo_reactive_osqp_solver_model_hpp__
#define __galileo_reactive_osqp_solver_model_hpp__

#include "galileo/reactive/qp-solver/qp-solver-base.hpp"
#include <OsqpEigen/OsqpEigen.h>

namespace galileo
{
    namespace reactive
    {
        template <typename Scalar, int Nx, int Ne, int Ni, int Options = 0>
        struct OSQPSolverModel : QPSolverModelBase<OSQPSolverModel<Scalar, Nx, Ne, Ni, Options>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using QPSolverData = QPSolverDataTpl<Scalar, Nx, Ne, Ni, Options>;

            void init(QPSolverData &&data)
            {
                solver_.settings()->setWarmStart(true);

                // Initialize the solver
                solver_.data()->setNumberOfVariables(data.NumDecisionVars());

                int numDecisionBounds = data.HasDecisionBounds() ? data.NumDecisionVars() : 0;

                // OSQP does not handle decision bounds or equality constraints separately, so we need to include them in the inequality matrix.
                solver_.data()->setNumberOfConstraints(data.NumEqualityConstraints() + data.NumInequalityConstraints() + numDecisionBounds);

                // Set the Hessian matrix
                solver_.data()->setHessianMatrix(data.H());

                // Set the linear term of the quadratic cost
                solver_.data()->setGradient(data.g());

                using NewIneq = typename math::VConMat<typename QPSolverData::MatExN, typename QPSolverData::MatIxN>::type;
                using NewIneqBound = typename math::VConVec<typename QPSolverData::VecNe, typename QPSolverData::VecNi>::type;

                // Convert the equality constraints to the inequality matrix
                NewIneq ineq = math::vertcat(std::move(data.A()), std::move(data.C()));
                NewIneqBound lb_in = math::vertcat(data.b(), std::move(data.l()));
                NewIneqBound ub_in = math::vertcat(data.b(), std::move(data.u()));

                // Set the linear constraints matrix
                if (data.HasDecisionBounds())
                {
                    using NewIneqWithDecBounds = typename math::VConMat<NewIneq, Eigen::MatrixXd>::type;
                    using NewIneqBoundWithDecBounds = typename math::VConVec<NewIneqBound, Eigen::MatrixXd>::type;

                    NewIneqWithDecBounds ineqWithDecBounds = math::vertcat(std::move(ineq), Eigen::MatrixXd::Identity(data.NumDecisionVars(), data.NumDecisionVars()));
                    NewIneqBoundWithDecBounds lbWithDecBounds = math::vertcat(std::move(lb_in), std::move(data.lb()));
                    NewIneqBoundWithDecBounds ubWithDecBounds = math::vertcat(std::move(ub_in), std::move(data.ub()));

                    solver_.data()->setLinearConstraintsMatrix(ineqWithDecBounds);
                    solver_.data()->setLowerBound(lbWithDecBounds);
                    solver_.data()->setUpperBound(ubWithDecBounds);
                }
                else
                {
                    solver_.data()->setLinearConstraintsMatrix(ineq);
                    solver_.data()->setLowerBound(lb_in);
                    solver_.data()->setUpperBound(ub_in);
                }

                // Initialize the solver
                solver_.initSolver();
            }

            void solve()
            {
                // Solve the QP
                solver_.solveProblem();
            }

            void update(QPSolverData &&data)
            {
                // Update the solver
                init(std::move(data));
            }

            template <typename VectorType>
            void getPrimal(Eigen::MatrixBase<VectorType> const &primal)
            {
                // Get the primal solution
                const_cast<Eigen::MatrixBase<VectorType> &>(primal) = solver_.getSolution();
            }

        private:
            OsqpEigen::Solver solver_;
        };

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_osqp_solver_model_hpp__