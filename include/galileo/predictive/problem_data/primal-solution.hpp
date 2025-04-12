#ifndef __galileo_predictive_problem_data_primal_solution_hpp__
#define __galileo_predictive_problem_data_primal_solution_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/references/mode-schedule-base.hpp"

#define GALILEO_PRIMAL_SOLUTION_BASIC_TYPEDEF(PrimalSolution)     \
    using VarScalar = typename traits<PrimalSolution>::VarScalar; \
    using NumScalar = typename traits<PrimalSolution>::NumScalar; \
    static constexpr int Options = traits<PrimalSolution>::Options;

#define GALILEO_PRIMAL_SOLUTION_TYPEDEF(PrimalSolution)                             \
    using ModeScheduleBase_t = typename traits<PrimalSolution>::ModeScheduleBase_t; \
    using VectorXv = typename traits<PrimalSolution>::VectorXv;                     \
    using VectorXn = typename traits<PrimalSolution>::VectorXn;

namespace galileo
{

    template <typename Derived>
    struct PrimalSolution : internal::CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PrimalSolutionDerived = typename traits<Derived>::PrimalSolutionDerived;
        GALILEO_PRIMAL_SOLUTION_BASIC_TYPEDEF(PrimalSolutionDerived);
        GALILEO_PRIMAL_SOLUTION_TYPEDEF(PrimalSolutionDerived);

        std::vector<NumScalar> time_trajectory_;
        std::vector<VectorXn> state_trajectory_;
        std::vector<VectorXn> input_trajectory_;
        std::vector<size_t> post_event_indices_;
        ModeScheduleBase_t mode_schedule_;

        // optimized controller

    }; // class PrimalSolution

} // namespace galileo

#endif // __galileo_predictive_problem_data_primal_solution_hpp__