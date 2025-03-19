#ifndef __galileo_predictive_mpc_base_hpp__
#define __galileo_predictive_mpc_base_hpp__

#include "galileo/predictive/fwd.hpp"

#include "galileo/predictive/trajectory.hpp"

#define GALILEO_MPC_BASIC_TYPEDEF(MPC)                 \
    using VarScalar = typename traits<MPC>::VarScalar; \
    using NumScalar = typename traits<MPC>::NumScalar; \
    static constexpr int Options = traits<MPC>::Options;

namespace galileo
{

    namespace predictive
    {

        template <class Derived>
        class MPCBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using MPCDerived = typename traits<Derived>::MPCDerived;
            GALILEO_MPC_BASIC_TYPEDEF(MPCDerived);

            void reset()
            {
                derived().reset();
            }

            template <typename StateVectorType>
            void run(NumScalar current_time, const Eigen::MatrixBase<StateVectorType> &x0)
            {
                derived().run(current_time, x0.derived());
            }

        }; // class MPCBase

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_mpc_base_hpp__