#ifndef __galileo_reactive_qp_solvers_qp_solver_data_base_hpp__
#define __galileo_reactive_qp_solvers_qp_solver_data_base_hpp__

#include "galileo/reactive/qp-solvers/qp-solver-base.hpp"
#include "galileo/math/concat.hpp"
#include <optional>

namespace galileo
{

    namespace reactive
    {

        // Unified QP solver data structure
        template <typename Scalar, int Nx, int Ne, int Ni, int Options = 0>
        class QPSolverDataTpl
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using MatNxN = Eigen::Matrix<Scalar, Nx, Nx, Options>;
            using VecNx = Eigen::Matrix<Scalar, Nx, 1, Options>;
            using MatExN = Eigen::Matrix<Scalar, Ne, Nx, Options>;
            using VecNe = Eigen::Matrix<Scalar, Ne, 1, Options>;
            using MatIxN = Eigen::Matrix<Scalar, Ni, Nx, Options>;
            using VecNi = Eigen::Matrix<Scalar, Ni, 1, Options>;

            QPSolverDataTpl() = default;

            // Perfect forwarding constructor
            template <typename H_in, typename g_in, typename A_in, typename b_in, typename C_in, typename l_in, typename u_in>
            QPSolverDataTpl(H_in &&H, g_in &&g, A_in &&A, b_in &&b, C_in &&C, l_in &&l, u_in &&u)
                : H_(math::toConcreteMatrix<H_in, MatNxN>(std::forward<H_in>(H))),
                  g_(math::toConcreteMatrix<g_in, VecNx>(std::forward<g_in>(g))),
                  A_(math::toConcreteMatrix<A_in, MatExN>(std::forward<A_in>(A))),
                  b_(math::toConcreteMatrix<b_in, VecNe>(std::forward<b_in>(b))),
                  C_(math::toConcreteMatrix<C_in, MatIxN>(std::forward<C_in>(C))),
                  l_(math::toConcreteMatrix<l_in, VecNi>(std::forward<l_in>(l))),
                  u_(math::toConcreteMatrix<u_in, VecNi>(std::forward<u_in>(u)))
            {
            }

            template <typename lb_in, typename ub_in>
            void SetDecisionBounds(const Eigen::MatrixBase<lb_in> &lb, const Eigen::MatrixBase<ub_in> &ub)
            {
                lb_ = lb;
                ub_ = ub;
            }

            int NumDecisionVars() const { return H_.rows(); }

            int NumEqualityConstraints() const { return A_.rows(); }

            int NumInequalityConstraints() const { return C_.rows(); }

            bool HasDecisionBounds() const { return lb_.has_value() && ub_.has_value(); }

            const MatNxN &H() const { return H_; }

            MatNxN &&H() { return std::move(H_); }

            const VecNx &g() const { return g_; }

            VecNx &&g() { return std::move(g_); }

            const MatExN &A() const { return A_; }

            MatExN &&A() { return std::move(A_); }

            const VecNe &b() const { return b_; }

            VecNe &&b() { return std::move(b_); }

            const MatIxN &C() const { return C_; }

            MatIxN &&C() { return std::move(C_); }

            const VecNi &l() const { return l_; }

            VecNi &&l() { return std::move(l_); }

            const VecNi &u() const { return u_; }

            VecNi &&u() { return std::move(u_); }

            const VecNx &lb() const { return *lb_; }

            VecNx &&lb() { return std::move(*lb_); }

            const VecNx &ub() const { return *ub_; }

            VecNx &&ub() { return std::move(*ub_); }

        protected:
            // Hessian and gradient
            MatNxN H_;
            VecNx g_;

            // Equality constraints: A x = b
            MatExN A_;
            VecNe b_;

            // Inequality constraints: l <= C x <= u
            MatIxN C_;
            VecNi l_;
            VecNi u_;

            // Optional bounds on x: lb <= x <= ub
            std::optional<VecNx> lb_;
            std::optional<VecNx> ub_;

        }; // class QPSolverDataTpl

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_qp_solvers_qp_solver_data_base_hpp__