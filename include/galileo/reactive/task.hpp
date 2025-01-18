#ifndef __galileo_reactive_task_hpp__
#define __galileo_reactive_task_hpp__

#include "galileo/reactive/fwd.hpp"
#include <type_traits>

namespace galileo
{
    namespace reactive
    {

        template <typename M1, typename M2>
        struct VConMat
        {
        private:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            static constexpr int R1 = Eigen::MatrixBase<M1>::RowsAtCompileTime;
            static constexpr int C1 = Eigen::MatrixBase<M1>::ColsAtCompileTime;
            static constexpr int R2 = Eigen::MatrixBase<M2>::RowsAtCompileTime;
            static constexpr int C2 = Eigen::MatrixBase<M2>::ColsAtCompileTime;

            // The new row count is either (R1+R2) if both are fixed, else Dynamic.
            static constexpr int Rows =
                (R1 != Eigen::Dynamic && R2 != Eigen::Dynamic)
                    ? (R1 + R2)
                    : Eigen::Dynamic;

            // If both have known, identical columns => keep that at compile time
            static constexpr bool sameFixedCols = (C1 != Eigen::Dynamic && C2 != Eigen::Dynamic && C1 == C2);
            static constexpr int Cols = sameFixedCols ? C1 : Eigen::Dynamic;

        public:
            // The underlying scalar type (e.g. double)
            using Scalar = typename M1::Scalar;

            // The resulting type
            using type = typename Eigen::Matrix<Scalar, Rows, Cols>;
        };

        template <typename V1, typename V2>
        struct VConVec
        {
        private:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            static constexpr int R1 = Eigen::MatrixBase<V1>::RowsAtCompileTime;
            static constexpr int R2 = Eigen::MatrixBase<V2>::RowsAtCompileTime;

            // The new row count is either (R1+R2) if both are fixed, else Dynamic.
            static constexpr int Rows =
                (R1 != Eigen::Dynamic && R2 != Eigen::Dynamic)
                    ? (R1 + R2)
                    : Eigen::Dynamic;

        public:
            // The underlying scalar type (e.g. double)
            using Scalar = typename V1::Scalar;

            // The resulting type
            using type = typename Eigen::Matrix<Scalar, Rows, 1>;
        };

        template <typename M1, typename M2>
        static typename VConMat<std::decay_t<M1>, std::decay_t<M2>>::type vertcat(M1 &&m1, M2 &&m2)
        {
            // 'std::decay_t<M1>' strips references and cv-qualifiers,
            // so if M1 is exactly some Eigen::Matrix<double,R,C> or expression,
            // we unify that type with VConMat.

            using M1Plain = std::decay_t<M1>;
            using M2Plain = std::decay_t<M2>;
            using ReturnType = typename VConMat<M1Plain, M2Plain>::type;

            // (1) Fully fixed in rows & cols
            if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                          ReturnType::ColsAtCompileTime != Eigen::Dynamic)
            {
                ReturnType res; // e.g. Matrix<double, R1+R2, C>
                res << m1, m2;  // one pass filling top/bottom
                return res;
            }
            // (2) Fixed cols, dynamic rows
            else if constexpr (ReturnType::RowsAtCompileTime == Eigen::Dynamic &&
                               ReturnType::ColsAtCompileTime != Eigen::Dynamic)
            {
                const int totalRows = m1.rows() + m2.rows();
                ReturnType res(totalRows, ReturnType::ColsAtCompileTime);
                res << m1, m2;
                return res;
            }
            // (3) Fixed rows, dynamic cols (unusual for vertical stacking, but included)
            else if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                               ReturnType::ColsAtCompileTime == Eigen::Dynamic)
            {
                const int totalCols = m1.cols(); // must match m2.cols() at runtime
                ReturnType res(ReturnType::RowsAtCompileTime, totalCols);
                res << m1, m2;
                return res;
            }
            // (4) Fully dynamic
            else
            {
                const int totalRows = m1.rows() + m2.rows();
                const int totalCols = m1.cols(); // must match m2.cols() at runtime
                ReturnType res(totalRows, totalCols);
                res << m1, m2;
                return res;
            }
        }

        namespace detail
        {

            // This helper decays an arbitrary Eigen expression or matrix (Derived)
            // into a concrete Matrix type MType.
            //
            // If Derived is exactly MType and an rvalue, we can move.
            // Otherwise, we construct (i.e. evaluate) a new MType.
            template <typename Derived, typename MType>
            EIGEN_STRONG_INLINE MType toConcreteMatrix(Derived &&x)
            {
                // Check if the decayed type of x is exactly MType:
                if constexpr (std::is_same_v<std::decay_t<Derived>, MType>)
                {
                    // We can move-construct MType if x is an rvalue of MType,
                    // or copy-construct if it's an lvalue. std::forward picks the right one.
                    return std::forward<Derived>(x);
                }
                else
                {
                    // Evaluate the expression or convert from a different type
                    // (e.g. from expression to Matrix, from block to Matrix, etc.)
                    return MType(std::forward<Derived>(x));
                }
            }

        } // end namespace detail

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
                : a_(detail::toConcreteMatrix<A_in, A_t>(std::forward<A_in>(a))),
                  b_(detail::toConcreteMatrix<B_in, B_t>(std::forward<B_in>(b))),
                  d_(detail::toConcreteMatrix<D_in, D_t>(std::forward<D_in>(d))),
                  f_(detail::toConcreteMatrix<F_in, F_t>(std::forward<F_in>(f)))
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
                a_.noalias() *= rhs;
                b_.noalias() *= rhs;
                d_.noalias() *= rhs;
                f_.noalias() *= rhs;

                return *this;
            }

            A_t a_;
            B_t b_;
            D_t d_;
            F_t f_;
        };

        template <typename T1, typename T2>
        static Task<typename VConMat<typename std::decay_t<T1>::A_t, typename std::decay_t<T2>::A_t>::type,
                    typename VConVec<typename std::decay_t<T1>::B_t, typename std::decay_t<T2>::B_t>::type,
                    typename VConMat<typename std::decay_t<T1>::D_t, typename std::decay_t<T2>::D_t>::type,
                    typename VConVec<typename std::decay_t<T1>::F_t, typename std::decay_t<T2>::F_t>::type>
        verticalAdd(T1 &&t1, T2 &&t2)
        {
            using T1dec = std::decay_t<T1>;
            using T2dec = std::decay_t<T1>;

            using AOut = typename VConMat<typename T1dec::A_t, typename T2dec::A_t>::type;
            using BOut = typename VConVec<typename T1dec::B_t, typename T2dec::B_t>::type;
            using DOut = typename VConMat<typename T1dec::D_t, typename T2dec::D_t>::type;
            using FOut = typename VConVec<typename T1dec::F_t, typename T2dec::F_t>::type;

            AOut a_new = vertcat(std::move(t1.a_), std::move(t2.a_));
            BOut b_new = vertcat(std::move(t1.b_), std::move(t2.b_));
            DOut d_new = vertcat(std::move(t1.d_), std::move(t2.d_));
            FOut f_new = vertcat(std::move(t1.f_), std::move(t2.f_));

            return Task<AOut, BOut, DOut, FOut>(a_new, b_new, d_new, f_new);
        }

    } // namespace reactive

} // namespace galileo

#endif // __galileo_reactive_task_hpp__