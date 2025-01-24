#ifndef __galileo_reactive_qp_solver_model_base_hpp__
#define __galileo_reactive_qp_solver_model_base_hpp__

#include "galileo/reactive/qp-solver/qp-solver-base.hpp"

namespace galileo
{
    namespace reactive
    {
        template <typename Derived>
        class QPSolverModelBase : CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            template <typename QPDataType>
            void init(QPDataType &&data)
            {
                derived().init(std::forward<QPDataType>(data));
            }

            void solve()
            {
                derived().solve();
            }

            template <typename QPDataType>
            void update(QPDataType &&data)
            {
                derived().update(std::forward<QPDataType>(data));
            }

            template <typename VectorType>
            void getPrimal(Eigen::MatrixBase<VectorType> const &primal)
            {
                derived().getPrimal(primal);
            }

        };

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_qp_solver_model_base_hpp__