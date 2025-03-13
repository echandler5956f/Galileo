#ifndef __galileo_predictive_solvers_solver_base_hpp__
#define __galileo_predictive_solvers_solver_base_hpp__

#include "galileo/predictive/solvers/fwd.hpp"
#include "galileo/predictive/optimal-control-problem.hpp"

#include <vector>
#include <memory>
#include <limits>

#define GALILEO_SOLVER_BASIC_TYPEDEF(Solver)              \
    using VarScalar = typename traits<Solver>::VarScalar; \
    using NumScalar = typename traits<Solver>::NumScalar; \
    static constexpr int Options = traits<Solver>::Options;

#define GALILEO_SOLVER_TYPEDEF(Solver)                                                \
    using OptimalControlProblem_t = typename traits<Solver>::OptimalControlProblem_t; \
    using VectorXvs = typename traits<Solver>::VectorXvs;                             \
    using VectorXns = typename traits<Solver>::VectorXns;                             \
    using Vector2ns = Eigen::Matrix<NumScalar, 2, 1, Options>;

namespace galileo
{

    namespace predictive
    {

        class CallbackAbstract; // forward declaration

        enum FeasibilityNorm
        {
            LInf = 0,
            L1
        };

        template <typename Derived>
        class SolverBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SolverDerived = typename traits<Derived>::SolverDerived;
            GALILEO_SOLVER_BASIC_TYPEDEF(SolverDerived);
            GALILEO_SOLVER_TYPEDEF(SolverDerived);

            bool solve(const std::vector<VectorXns> &init_xs,
                       const std::vector<VectorXns> &init_us,
                       const std::size_t max_iter = 100, const bool is_feasible = false,
                       const NumScalar reg_init = NAN)
            {
                return derived().solve(init_xs, init_us, max_iter, is_feasible, reg_init);
            }

            void computeDirection(const bool recalc)
            {
                derived().computeDirection(recalc);
            }

            NumScalar tryStep(const NumScalar steplength)
            {
                return derived().tryStep(steplength);
            }

            NumScalar stoppingCriteria()
            {
                return derived().stoppingCriteria();
            }

            Vector2ns expectedImprovement()
            {
                return derived().expectedImprovement();
            }

            void resizeData()
            {
                derived().resizeData();
            }

            NumScalar computeDynamicFeasibility();

            NumScalar computeEqualityFeasibility();

            NumScalar computeInequalityFeasibility();

            void setCandidate(const std::vector<VectorXns> &xs_warm,
                              const std::vector<VectorXns> &us_warm,
                              bool is_feasible = false);

            void set_xs(const std::vector<VectorXns> &xs);

            void set_us(const std::vector<VectorXns> &us);

        protected:
            inline SolverBase()
            {
            }

            inline SolverBase(const SolverBase &clone)
            {
                *this = clone;
            }

            inline SolverBase &operator=(const SolverBase &clone)
            {
                return *this;
            }

            OptimalControlProblem_t ocp_;

            std::vector<VectorXns> x_;
            std::vector<VectorXns> u_;

            std::vector<VectorXns> fs_;

            std::vector<std::shared_ptr<CallbackAbstract>> callbacks_;
            bool is_feasible_;
            bool was_feasible_;

            NumScalar cost_;
            NumScalar merit_;
            NumScalar stop_;
            Vector2ns d_;
            NumScalar dV_;
            NumScalar dPhi_;
            NumScalar dVexp_;
            NumScalar dPhiexp_;
            NumScalar dfeas_;
            NumScalar feas_;

            NumScalar ffeas_;
            NumScalar gfeas_;

            NumScalar hfeas_;

            NumScalar ffeas_try_;

            NumScalar gfeas_try_;

            NumScalar hfeas_try_;

            NumScalar preg_;
            NumScalar dreg_;

            NumScalar steplength_;
            NumScalar th_acceptstep_;
            NumScalar th_stop_;
            NumScalar th_gaptol_;
            enum FeasibilityNorm feasnorm_;

            std::size_t iter_;
            NumScalar tmp_feas_;
            std::vector<VectorXns> g_adj_;

        }; // class SolverBase

    } // namespace predictive

} // namespace galileo

#include "galileo/predictive/solvers/solver-base.hxx"

#endif // __galileo_predictive_solvers_solver_base_hpp__