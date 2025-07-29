#ifndef __galileo_solvers_solver_base_hpp__
#define __galileo_solvers_solver_base_hpp__

#include "galileo/solvers/fwd.hpp"

namespace galileo
{

    template <typename Derived>
    class SolverBase
        : public internal::CRTP<Derived>
    {
    public:
        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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

    }; // class SolverBase

} // namespace galileo

#endif // __galileo_solvers_solver_base_hpp__
