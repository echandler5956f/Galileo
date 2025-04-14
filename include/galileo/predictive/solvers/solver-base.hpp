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
    static constexpr FeasibilityNormOptions feasnorm_ = traits<Solver>::feasnorm_;                \
    using VectorXv = typename traits<Solver>::VectorXv;                               \
    using VectorXn = typename traits<Solver>::VectorXn;

namespace galileo
{

    enum class FeasibilityNormOptions
    {
        LInf = 0,
        L1 = 1
    }; // enum FeasibilityNormOptions

    template <typename Derived>
    class SolverBase : internal::CRTP<SolverBase<Derived>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using SolverDerived = typename traits<Derived>::SolverDerived;
        GALILEO_SOLVER_BASIC_TYPEDEF(SolverDerived);
        GALILEO_SOLVER_TYPEDEF(SolverDerived);

        bool solve(const std::vector<VectorXn> &init_xs,
                   const std::vector<VectorXn> &init_us,
                   const std::size_t maxiter = 100,
                   const bool is_feasible = false,
                   const NumScalar init_reg = NAN)
        {
            return this->derived().solve(init_xs, init_us, maxiter, is_feasible, init_reg);
        }

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

        void resizeData();

        NumScalar computeDynamicFeasibility();

        NumScalar computeInequalityFeasibility();

        NumScalar computeEqualityFeasibility();

        void setCandidate(const std::vector<VectorXn> &xs_warm, const std::vector<VectorXn> &us_warm, const bool is_feasible = false);

        OptimalControlProblem_t ocp_;

        std::vector<VectorXn> xs_;
        std::vector<VectorXn> us_;
        std::vector<VectorXn> fs_;

        bool is_feasible_;
        bool was_feasible_;
        NumScalar cost_;

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

        std::size_t iter_;
        NumScalar tmp_feas_;
        std::vector<VectorXns> g_adj_;

    }; // class SolverBase

} // namespace galileo

// #include "galileo/predictive/solvers/solver-base.hxx"

#endif // __galileo_predictive_solvers_solver_base_hpp__