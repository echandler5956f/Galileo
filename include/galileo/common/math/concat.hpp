#ifndef __galileo_common_linalg_concat_hpp__
#define __galileo_common_linalg_concat_hpp__

#include "galileo/common/fwd.hpp"

namespace galileo
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

        // Ensure the number of columns agree at runtime.
        GALILEO_ASSERT(m1.cols() == m2.cols(), "VConMat: Matrices must have the same number of columns.");

        // Grab the runtime number of rows.
        const int rows1 = m1.rows();
        const int rows2 = m2.rows();

        // (A) If one (or both) matrices is empty, simply return the nonempty one.
        if (rows1 == 0 && rows2 == 0)
        {
            // Both are empty: return an empty matrix with the proper number of columns.
            return ReturnType(0, m1.cols());
        }
        else if (rows1 == 0)
        {
            // m1 is empty: return m2 converted to the ReturnType.
            return ReturnType(m2);
        }
        else if (rows2 == 0)
        {
            // m2 is empty: return m1 converted to the ReturnType.
            return ReturnType(m1);
        }

        // (B) Neither matrix is empty. Now do the normal vertical concatenation.
        if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                      ReturnType::ColsAtCompileTime != Eigen::Dynamic)
        {
            // Fully fixed-size: let Eigen's comma initializer fill the result.
            ReturnType res; // (size is fixed at compile time)
            res << m1, m2;  // fills res in one pass (top rows from m1, bottom from m2)
            return res;
        }
        else if constexpr (ReturnType::RowsAtCompileTime == Eigen::Dynamic &&
                           ReturnType::ColsAtCompileTime != Eigen::Dynamic)
        {
            const int totalRows = rows1 + rows2;
            ReturnType res(totalRows, ReturnType::ColsAtCompileTime);
            res << m1, m2;
            return res;
        }
        else if constexpr (ReturnType::RowsAtCompileTime != Eigen::Dynamic &&
                           ReturnType::ColsAtCompileTime == Eigen::Dynamic)
        {
            // Unusual for vertical stacking (fixed rows, dynamic columns).
            const int totalCols = m1.cols(); // must match m2.cols() at runtime
            ReturnType res(ReturnType::RowsAtCompileTime, totalCols);
            res << m1, m2;
            return res;
        }
        else // (Fully dynamic)
        {
            const int totalRows = rows1 + rows2;
            const int totalCols = m1.cols(); // again, m1.cols() == m2.cols()
            ReturnType res(totalRows, totalCols);
            res << m1, m2;
            return res;
        }
    }

} // namespace galileo

#endif // __galileo_common_linalg_concat_hpp__
