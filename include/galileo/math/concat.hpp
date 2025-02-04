#ifndef __galileo_math_concat_hpp__
#define __galileo_math_concat_hpp__

#include "galileo/math/fwd.hpp"

namespace galileo
{
    namespace math
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

        }; // struct VConVec

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

    } // namespace math

} // namespace galileo

#endif // __galileo_math_concat_hpp__