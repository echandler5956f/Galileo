#ifndef __galileo_reactive_task_hpp__
#define __galileo_reactive_task_hpp__

#include "galileo/reactive/fwd.hpp"
#include "galileo/math/matrix.hpp"
#include "galileo/math/concat.hpp"
#include <type_traits>

namespace galileo
{
    namespace reactive
    {

        template <typename AType, typename BType, typename DType, typename FType>
        class Task
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using A_t = AType;
            using B_t = BType;
            using D_t = DType;
            using F_t = FType;

            using Scalar = typename A_t::Scalar;

            // Perfect-forwarding constructor
            template <typename A_in, typename B_in, typename D_in, typename F_in>
            Task(A_in &&a, B_in &&b, D_in &&d, F_in &&f)
                : a_(math::toConcreteMatrix<A_in, A_t>(std::forward<A_in>(a))),
                  b_(math::toConcreteMatrix<B_in, B_t>(std::forward<B_in>(b))),
                  d_(math::toConcreteMatrix<D_in, D_t>(std::forward<D_in>(d))),
                  f_(math::toConcreteMatrix<F_in, F_t>(std::forward<F_in>(f)))
            {
            }

            explicit Task(const size_t numDecisionVars) : Task(A_t::Zero(numDecisionVars, numDecisionVars),
                                                               B_t::Zero(numDecisionVars, 1),
                                                               D_t::Zero(numDecisionVars, numDecisionVars),
                                                               F_t::Zero(numDecisionVars, 1))
            {
            }

            // Provide an operator* that scales the A, B, D, F matrices by a scalar.
            Task<A_t, B_t, D_t, F_t> &operator*(const Scalar rhs)
            {
                a_.noalias() = a_ * rhs;
                b_.noalias() = b_ * rhs;
                d_.noalias() = d_ * rhs;
                f_.noalias() = f_ * rhs;

                return *this;
            }

            A_t a_;
            B_t b_;
            D_t d_;
            F_t f_;
        };

        template <typename Scalar>
        using TaskDefault = Task<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>, Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>;

        template <typename T1, typename T2>
        static Task<typename math::VConMat<typename std::decay_t<T1>::A_t, typename std::decay_t<T2>::A_t>::type,
                    typename math::VConVec<typename std::decay_t<T1>::B_t, typename std::decay_t<T2>::B_t>::type,
                    typename math::VConMat<typename std::decay_t<T1>::D_t, typename std::decay_t<T2>::D_t>::type,
                    typename math::VConVec<typename std::decay_t<T1>::F_t, typename std::decay_t<T2>::F_t>::type>
        taskVertcat(T1 &&t1, T2 &&t2)
        {
            using T1dec = std::decay_t<T1>;
            using T2dec = std::decay_t<T2>;

            using AOut = typename math::VConMat<typename T1dec::A_t, typename T2dec::A_t>::type;
            using BOut = typename math::VConVec<typename T1dec::B_t, typename T2dec::B_t>::type;
            using DOut = typename math::VConMat<typename T1dec::D_t, typename T2dec::D_t>::type;
            using FOut = typename math::VConVec<typename T1dec::F_t, typename T2dec::F_t>::type;

            AOut a_new = math::vertcat(std::move(t1.a_), std::move(t2.a_));
            BOut b_new = math::vertcat(std::move(t1.b_), std::move(t2.b_));
            DOut d_new = math::vertcat(std::move(t1.d_), std::move(t2.d_));
            FOut f_new = math::vertcat(std::move(t1.f_), std::move(t2.f_));
            return Task<AOut, BOut, DOut, FOut>(a_new, b_new, d_new, f_new);
        }

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_task_hpp__