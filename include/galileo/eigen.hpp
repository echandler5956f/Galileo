#ifndef __galileo_eigen_hpp__
#define __galileo_eigen_hpp__

#include "galileo/dimensions.hpp"
#include <Eigen/Core>
#include <concepts>

namespace galileo
{

    // A concept to constrain templates to Eigen types
    template <typename Derived>
    concept EigenRowVector = (Derived::RowsAtCompileTime == 1 && Derived::ColsAtCompileTime != 1);

    template <typename Derived>
    concept EigenColVector = (Derived::ColsAtCompileTime == 1 && Derived::RowsAtCompileTime != 1);

    template <typename Derived>
    concept EigenVector = EigenRowVector<Derived> || EigenColVector<Derived>;

    template <typename Derived>
    concept EigenMatrix = !EigenVector<Derived>;

    template <int R = Eigen::Dynamic, int C = Eigen::Dynamic>
    struct AccessDispatcher
    {
        static_assert((R == Eigen::Dynamic || R > 0) &&
                          (C == Eigen::Dynamic || C > 0),
                      "Compile-time sizes must be positive or Eigen::Dynamic");

        using RowDim = Dimension<R>;
        using ColDim = Dimension<C>;

        // Vectors (only enabled when expression is a vector)
        static constexpr int LenCT = (R == 1 ? C : (C == 1 ? R : -1));

        template <typename Derived>
            requires EigenVector<Derived>
        static auto segment(Eigen::MatrixBase<Derived> &vec,
                            Eigen::Index start,
                            const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template segment<LenCT>(start);
            else
                return vec.segment(start, len);
        }

        template <typename Derived>
            requires EigenVector<Derived>
        static auto segment(const Eigen::MatrixBase<Derived> &vec,
                            Eigen::Index start,
                            const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template segment<LenCT>(start);
            else
                return vec.segment(start, len);
        }

        template <typename Derived>
            requires EigenVector<Derived>
        static auto head(Eigen::MatrixBase<Derived> &vec,
                         const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template head<LenCT>();
            else
                return vec.head(len);
        }

        template <typename Derived>
            requires EigenVector<Derived>
        static auto head(const Eigen::MatrixBase<Derived> &vec,
                         const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template head<LenCT>();
            else
                return vec.head(len);
        }

        template <typename Derived>
            requires EigenVector<Derived>
        static auto tail(Eigen::MatrixBase<Derived> &vec,
                         const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template tail<LenCT>();
            else
                return vec.tail(len);
        }

        template <typename Derived>
            requires EigenVector<Derived>
        static auto tail(const Eigen::MatrixBase<Derived> &vec,
                         const Dimension<LenCT> &len = Dimension<LenCT>())
        {
            if constexpr (LenCT != Eigen::Dynamic)
                return vec.template tail<LenCT>();
            else
                return vec.tail(len);
        }

        // Matrices (rectangular or square)
        template <typename Derived>
        static auto block(Eigen::MatrixBase<Derived> &mat,
                          Eigen::Index row0,
                          Eigen::Index col0,
                          const RowDim &rsz = RowDim(),
                          const ColDim &csz = ColDim())
        {
            constexpr bool RowFixed = (R != Eigen::Dynamic);
            constexpr bool ColFixed = (C != Eigen::Dynamic);

            if constexpr (RowFixed && ColFixed)
                return mat.template block<R, C>(row0, col0);
            else if constexpr (RowFixed && !ColFixed)
                return mat.template block<R, Eigen::Dynamic>(row0, col0, R, csz);
            else if constexpr (!RowFixed && ColFixed)
                return mat.template block<Eigen::Dynamic, C>(row0, col0, rsz, C);
            else
                return mat.block(row0, col0, rsz, csz);
        }

        template <typename Derived>
        static auto block(const Eigen::MatrixBase<Derived> &mat,
                          Eigen::Index row0,
                          Eigen::Index col0,
                          const RowDim &rsz = RowDim(),
                          const ColDim &csz = ColDim())
        {
            constexpr bool RowFixed = (R != Eigen::Dynamic);
            constexpr bool ColFixed = (C != Eigen::Dynamic);

            if constexpr (RowFixed && ColFixed)
                return mat.template block<R, C>(row0, col0);
            else if constexpr (RowFixed && !ColFixed)
                return mat.template block<R, Eigen::Dynamic>(row0, col0, R, csz);
            else if constexpr (!RowFixed && ColFixed)
                return mat.template block<Eigen::Dynamic, C>(row0, col0, rsz, C);
            else
                return mat.block(row0, col0, rsz, csz);
        }

        template <typename Mat>
        using BlockXpr = Eigen::Block<Mat, R, C>;
    }; // struct AccessDispatcher

    template <int N = Eigen::Dynamic>
    using VectorDispatcher = AccessDispatcher<N, 1>;

    template <int R = Eigen::Dynamic, int C = Eigen::Dynamic>
    using MatrixDispatcher = AccessDispatcher<R, C>;

    // Helper alias that automatically corrects the storage order for 1xN or Nx1
    // This is necessary because Eigen column vectors must be stored in column-major order,
    // while row vectors must be stored in row-major order. If a matrix is fed two constants,
    // where either the row or column dimensions may or may not have a size of 1,
    // it is convenient to have the storage order be automatically determined rather than
    // requiring the user to specify it.
    template <typename Scalar, int Rows, int Cols, int DefaultOptions>
    struct MatrixTpl
    {
        using nominal_type = Eigen::Matrix<Scalar, Rows, Cols, DefaultOptions>;

        static constexpr bool isRowVector = EigenRowVector<nominal_type>;
        static constexpr bool isColVector = EigenColVector<nominal_type>;

        static constexpr int OptionsFixed =
            isRowVector   ? ((DefaultOptions & ~Eigen::RowMajor) | Eigen::RowMajor)
            : isColVector ? (DefaultOptions & ~Eigen::RowMajor)
                          : DefaultOptions;

        static_assert(!(isRowVector && !(OptionsFixed & Eigen::RowMajor)),
                      "A 1xN Eigen matrix must be row-major.");
        static_assert(!(isColVector && (OptionsFixed & Eigen::RowMajor)),
                      "An Nx1 Eigen matrix must be column-major.");

        using type = Eigen::Matrix<Scalar, Rows, Cols, OptionsFixed>;
    }; // struct MatrixTpl

    template <typename Scalar, int Rows, int Cols, int DefaultOptions = Eigen::ColMajor>
    using Matrix = typename MatrixTpl<Scalar, Rows, Cols, DefaultOptions>::type;

} // namespace galileo

#endif // __galileo_eigen_hpp__
