#ifndef __galileo_predictive_solvers_solver_base_hxx__
#define __galileo_predictive_solvers_solver_base_hxx__

// #include "galileo/predictive/solvers/solver-base.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename Derived>
        SolverBase<Derived>::NumScalar SolverBase<Derived>::computeDynamicFeasibility()
        {
        }

        template <typename Derived>
        SolverBase<Derived>::NumScalar SolverBase<Derived>::computeEqualityFeasibility()
        {
        }

        template <typename Derived>
        SolverBase<Derived>::NumScalar SolverBase<Derived>::computeInequalityFeasibility()
        {
        }

        template <typename Derived>
        void SolverBase<Derived>::setCandidate(const std::vector<SolverBase<Derived>::VectorXns> &xs_warm,
                                               const std::vector<SolverBase<Derived>::VectorXns> &us_warm,
                                               bool is_feasible = false)
        {
        }

        template <typename Derived>
        void SolverBase<Derived>::set_xs(const std::vector<SolverBase<Derived>::VectorXns> &xs)
        {
        }

        template <typename Derived>
        void SolverBase<Derived>::set_us(const std::vector<SolverBase<Derived>::VectorXns> &us)
        {
        }

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_solvers_solver_base_hxx__