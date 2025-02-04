#ifndef __galileo_reactive_qp_solvers_proxqp_solver_hpp__
#define __galileo_reactive_qp_solvers_proxqp_solver_hpp__

#include "galileo/reactive/qp-solvers/qp-solver-base.hpp"
#include <proxsuite/proxqp/dense/dense.hpp>

namespace galileo
{

    namespace reactive
    {

        template <typename Scalar, int Nx, int Ne, int Ni, int Options = 0>
        class ProxQPSolver : QPSolverModelBase<ProxQPSolver<Scalar, Nx, Ne, Ni, Options>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using QPSolverData = QPSolverDataTpl<Scalar, Nx, Ne, Ni, Options>;

            void init(QPSolverData &&data)
            {
                // Initialize the solver
                if (data.HasDecisionBounds())
                {
                    using NewIneq = typename math::VConMat<typename QPSolverData::MatIxN, Eigen::MatrixXd>::type;
                    using NewIneqBound = typename math::VConVec<typename QPSolverData::VecNi, Eigen::MatrixXd>::type;

                    NewIneq ineq = math::vertcat(std::move(data.C()), Eigen::MatrixXd::Identity(data.NumDecisionVars(), data.NumDecisionVars()));
                    NewIneqBound lb_in = math::vertcat(std::move(data.l()), std::move(data.lb()));
                    NewIneqBound ub_in = math::vertcat(std::move(data.u()), std::move(data.ub()));

                    solver_.init(data.H, data.g, data.A, data.b, ineq, lb_in, ub_in);
                }
                else
                {
                    solver_.init(data.H, data.g, data.A, data.b, data.C, data.l, data.u);
                }
            }

            void solve()
            {
                // Solve the QP
                solver_.solve();
            }

            void update(QPSolverData &&data)
            {
                solver_.settings.initial_guess = proxsuite::proxqp::InitialGuessStatus::WARM_START_WITH_PREVIOUS_RESULT;

                // Update the solver
                if (data.HasDecisionBounds())
                {
                    using NewIneq = typename math::VConMat<typename QPSolverData::MatIxN, Eigen::MatrixXd>::type;
                    using NewIneqBound = typename math::VConVec<typename QPSolverData::VecNi, Eigen::MatrixXd>::type;

                    NewIneq ineq = math::vertcat(std::move(data.C()), Eigen::MatrixXd::Identity(data.NumDecisionVars(), data.NumDecisionVars()));
                    NewIneqBound lb_in = math::vertcat(std::move(data.l()), std::move(data.lb()));
                    NewIneqBound ub_in = math::vertcat(std::move(data.u()), std::move(data.ub()));

                    solver_.update(data.H, data.g, data.A, data.b, ineq, lb_in, ub_in);
                }
                else
                {
                    solver_.update(data.H, data.g, data.A, data.b, data.C, data.l, data.u);
                }
            }

            template <typename VectorType>
            void getPrimal(Eigen::MatrixBase<VectorType> const &primal)
            {
                // Get the primal solution
                const_cast<Eigen::MatrixBase<VectorType> &>(primal) = solver_.results.x;
            }

        private:
            proxsuite::proxqp::dense::QP<double> solver_;

        }; // class ProxQPSolverModel

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_qp_solvers_proxqp_solver_hpp__