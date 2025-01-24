#ifndef __galileo_reactive_qpoases_solver_model_hpp__
#define __galileo_reactive_qpoases_solver_model_hpp__

#include "galileo/reactive/qp-solver/qp-solver-base.hpp"
#include <qpOASES.hpp>

namespace galileo
{
    namespace reactive
    {
        template <typename Scalar, int Nx, int Ne, int Ni, int Options = 0>
        struct QPOASESSolverModel : QPSolverModelBase<QPOASESSolverModel<Scalar, Nx, Ne, Ni, Options>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using QPSolverData = QPSolverDataTpl<Scalar, Nx, Ne, Ni, Options>;

            void init(QPSolverData &&data)
            {
                qpOASES::QProblem solver = qpOASES::QProblem(data.NumDecisionVars(), data.NumEqualityConstraints() + data.NumInequalityConstraints());

                qpOASES::Options options;
                options.setToMPC();
                options.printLevel = qpOASES::PL_LOW;
                options.enableEqualities = qpOASES::BT_TRUE;
                qpProblem.setOptions(options);
                int nWsr = 20;

                using NewIneq = typename math::VConMat<typename QPSolverData::MatExN, typename QPSolverData::MatIxN>::type;
                using NewIneqBound = typename math::VConVec<typename QPSolverData::VecNe, typename QPSolverData::VecNi>::type;

                // Convert the equality constraints to the inequality matrix
                NewIneq ineq = math::vertcat(std::move(data.A()), std::move(data.C()));
                NewIneqBound lb_in = math::vertcat(data.b(), std::move(data.l()));
                NewIneqBound ub_in = math::vertcat(data.b(), std::move(data.u()));

                if (data.HasDecisionBounds())
                    solver.init(data.H().data(), data.g().data(), ineq.data(), data.lb().data(), data.ub().data(), lb_in.data(), ub_in.data(), nWsr);
                else
                    solver.init(data.H().data(), data.g().data(), ineq.data(), nullptr, nullptr, lb_in.data(), ub_in.data(), nWsr);

                solver.getPrimalSolution(primal_.data());
            }

            void solve()
            {
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
                const_cast<Eigen::MatrixBase<VectorType> &>(primal) = primal_;
            }

        private:
            Eigen::Matrix<Scalar, Nx, 1, Options> primal_;
        };

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_qpoases_solver_model_hpp__