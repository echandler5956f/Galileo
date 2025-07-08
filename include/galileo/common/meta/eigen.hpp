#ifndef __galileo_common_meta_eigen_hpp__
#define __galileo_common_meta_eigen_hpp__

#include <cassert>
#include <cmath>
#include <type_traits>

#include <Eigen/Core>

#include "galileo/common/meta/concepts.hpp"
#include "galileo/common/meta/dimensions.hpp"

namespace galileo
{

    namespace detail
    {
        // Helper alias that automatically corrects the storage order for 1xN or Nx1
        // This is necessary because Eigen column vectors must be stored in column-major order,
        // while row vectors must be stored in row-major order. If a matrix is fed two constants,
        // where either the row or column dimensions may or may not have a size of 1,
        // it is convenient to have the storage order be automatically determined rather than
        // requiring the user to specify it.
        template <typename Scalar, int Rows, int Cols, int DefaultOptions, int MaxRows = Rows, int MaxCols = Cols>
        struct MatrixTpl
        {
            using NominalType = Eigen::Matrix<Scalar, Rows, Cols, DefaultOptions, MaxRows, MaxCols>;

            static constexpr bool IsRowVector = IsEigenRowVector<NominalType>;
            static constexpr bool IsColVector = IsEigenColVector<NominalType>;

            static constexpr int OptionsFixed =
                IsRowVector   ? ((DefaultOptions & ~Eigen::RowMajor) | Eigen::RowMajor)
                : IsColVector ? (DefaultOptions & ~Eigen::RowMajor)
                              : DefaultOptions;

            static_assert(!(IsRowVector && !(OptionsFixed & Eigen::RowMajor)),
                          "A 1xN Eigen matrix must be row-major.");
            static_assert(!(IsColVector && (OptionsFixed & Eigen::RowMajor)),
                          "An Nx1 Eigen matrix must be column-major.");

            using Type = Eigen::Matrix<Scalar, Rows, Cols, OptionsFixed, MaxRows, MaxCols>;
        }; // struct MatrixTpl

    } // namespace detail

    // All matrices in Galileo are declared with galileo::Matrix (also available as Eigen::GMatrix).
    // It is an Eigen::Matrix under the hood; we only correct the storage order automatically when the
    // static shape degenerates into a row- or column-vector. You get the full Eigen API and the right
    // storage order without thinking about it.
    template <typename Scalar, int Rows, int Cols, int DefaultOptions = Eigen::ColMajor, int MaxRows = Rows, int MaxCols = Cols>
    using Matrix = typename detail::MatrixTpl<Scalar, Rows, Cols, DefaultOptions, MaxRows, MaxCols>::Type;

} // namespace galileo

namespace Eigen
{
    // Eigen::GMatrix is just an alias that picks a safe storage-order when R==1 or C==1;
    // otherwise it behaves exactly like Eigen::Matrix.
    template <typename Scalar, int Rows, int Cols, int DefaultOptions = Eigen::ColMajor, int MaxRows = Rows, int MaxCols = Cols>
    using GMatrix = galileo::Matrix<Scalar, Rows, Cols, DefaultOptions, MaxRows, MaxCols>;

} // namespace Eigen

namespace galileo
{

    // Segment accessors
    template <typename VecType, typename LenType>
    static auto segmentImpl(VecType &&vec, Eigen::Index start, const LenType &len)
    {
        if constexpr (std::is_integral_v<LenType>)
        {
            return vec.segment(start, len);
        }
        else
        {
            // DimensionTpl case
            if constexpr (LenType::IsDynamic)
            {
                return vec.segment(start, len.value());
            }
            else
            {
                return vec.template segment<LenType::Value>(start);
            }
        }
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto segment(Eigen::MatrixBase<Derived> &vec,
                        Eigen::Index start)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template segment<compile_time_len>(start);
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto segment(const Eigen::MatrixBase<Derived> &vec,
                        Eigen::Index start)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template segment<compile_time_len>(start);
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto segment(Eigen::MatrixBase<Derived> &vec,
                        Eigen::Index start,
                        const LenType &len)
    {
        return segmentImpl(vec, start, len);
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto segment(const Eigen::MatrixBase<Derived> &vec,
                        Eigen::Index start,
                        const LenType &len)
    {
        return segmentImpl(vec, start, len);
    }

    // Head accessors
    template <typename VecType, typename LenType>
    static auto headImpl(VecType &&vec, const LenType &len)
    {
        if constexpr (std::is_integral_v<LenType>)
        {
            return vec.head(len);
        }
        else
        {
            if constexpr (LenType::IsDynamic)
            {
                return vec.head(len.value());
            }
            else
            {
                return vec.template head<LenType::Value>();
            }
        }
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto head(Eigen::MatrixBase<Derived> &vec)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template head<compile_time_len>();
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto head(const Eigen::MatrixBase<Derived> &vec)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template head<compile_time_len>();
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto head(Eigen::MatrixBase<Derived> &vec, const LenType &len)
    {
        return headImpl(vec, len);
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto head(const Eigen::MatrixBase<Derived> &vec, const LenType &len)
    {
        return headImpl(vec, len);
    }

    // Tail accessors
    template <typename VecType, typename LenType>
    static auto tailImpl(VecType &&vec, const LenType &len)
    {
        if constexpr (std::is_integral_v<LenType>)
        {
            return vec.tail(len);
        }
        else
        {
            if constexpr (LenType::IsDynamic)
            {
                return vec.tail(len.value());
            }
            else
            {
                return vec.template tail<LenType::Value>();
            }
        }
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto tail(Eigen::MatrixBase<Derived> &vec)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template tail<compile_time_len>();
    }

    template <auto Len, typename Derived>
        requires IsEigenVector<Derived>
    static auto tail(const Eigen::MatrixBase<Derived> &vec)
    {
        constexpr int compile_time_len = extract_compile_time_value<Len>::Value;
        return vec.template tail<compile_time_len>();
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto tail(Eigen::MatrixBase<Derived> &vec, const LenType &len)
    {
        return tailImpl(vec, len);
    }

    template <typename Derived, typename LenType>
        requires IsEigenVector<Derived>
    static auto tail(const Eigen::MatrixBase<Derived> &vec, const LenType &len)
    {
        return tailImpl(vec, len);
    }

    // Block accessors
    template <typename MatType, typename RowSizeType, typename ColSizeType>
    static auto blockImpl(MatType &&mat, int start_row, int start_col,
                          const RowSizeType &row_size, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            return mat.block(start_row, start_col, row_size, col_size);
        }
        else if constexpr (std::is_integral_v<RowSizeType> && !std::is_integral_v<ColSizeType>)
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.block(start_row, start_col, row_size, col_size.value());
            }
            else
            {
                return mat.template block<Eigen::Dynamic, ColSizeType::Value>(start_row, start_col, row_size, ColSizeType::Value);
            }
        }
        else if constexpr (!std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.block(start_row, start_col, row_size.value(), col_size);
            }
            else
            {
                return mat.template block<RowSizeType::Value, Eigen::Dynamic>(start_row, start_col, RowSizeType::Value, col_size);
            }
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic && ColSizeType::IsDynamic)
            {
                return mat.block(start_row, start_col, row_size.value(), col_size.value());
            }
            else if constexpr (RowSizeType::IsDynamic && ColSizeType::IsFixed)
            {
                return mat.template block<Eigen::Dynamic, ColSizeType::Value>(start_row, start_col, row_size.value(), ColSizeType::Value);
            }
            else if constexpr (RowSizeType::IsFixed && ColSizeType::IsDynamic)
            {
                return mat.template block<RowSizeType::Value, Eigen::Dynamic>(start_row, start_col, RowSizeType::Value, col_size.value());
            }
            else
            {
                return mat.template block<RowSizeType::Value, ColSizeType::Value>(start_row, start_col);
            }
        }
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto block(Eigen::MatrixBase<Derived> &mat, int start_row, int start_col)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template block<compile_time_rows, compile_time_cols>(start_row, start_col);
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto block(const Eigen::MatrixBase<Derived> &mat, int start_row, int start_col)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template block<compile_time_rows, compile_time_cols>(start_row, start_col);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto block(Eigen::MatrixBase<Derived> &mat, int start_row, int start_col,
                      const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return blockImpl(mat, start_row, start_col, row_size, col_size);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto block(const Eigen::MatrixBase<Derived> &mat, int start_row, int start_col,
                      const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return blockImpl(mat, start_row, start_col, row_size, col_size);
    }

    // Top-left corner accessors
    template <typename MatType, typename RowSizeType, typename ColSizeType>
    static auto topLeftCornerImpl(MatType &&mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            return mat.topLeftCorner(row_size, col_size);
        }
        else if constexpr (std::is_integral_v<RowSizeType> && !std::is_integral_v<ColSizeType>)
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.topLeftCorner(row_size, col_size.value());
            }
            else
            {
                return mat.template topLeftCorner<Eigen::Dynamic, ColSizeType::Value>(row_size, ColSizeType::Value);
            }
        }
        else if constexpr (!std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.topLeftCorner(row_size.value(), col_size);
            }
            else
            {
                return mat.template topLeftCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size);
            }
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic && ColSizeType::IsDynamic)
            {
                return mat.topLeftCorner(row_size.value(), col_size.value());
            }
            else if constexpr (RowSizeType::IsDynamic && ColSizeType::IsFixed)
            {
                return mat.template topLeftCorner<Eigen::Dynamic, ColSizeType::Value>(row_size.value(), ColSizeType::Value);
            }
            else if constexpr (RowSizeType::IsFixed && ColSizeType::IsDynamic)
            {
                return mat.template topLeftCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size.value());
            }
            else
            {
                return mat.template topLeftCorner<RowSizeType::Value, ColSizeType::Value>();
            }
        }
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto topLeftCorner(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template topLeftCorner<compile_time_rows, compile_time_cols>();
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto topLeftCorner(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template topLeftCorner<compile_time_rows, compile_time_cols>();
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto topLeftCorner(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return topLeftCornerImpl(mat, row_size, col_size);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto topLeftCorner(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return topLeftCornerImpl(mat, row_size, col_size);
    }

    // Top-right corner accessors
    template <typename MatType, typename RowSizeType, typename ColSizeType>
    static auto topRightCornerImpl(MatType &&mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            return mat.topRightCorner(row_size, col_size);
        }
        else if constexpr (std::is_integral_v<RowSizeType> && !std::is_integral_v<ColSizeType>)
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.topRightCorner(row_size, col_size.value());
            }
            else
            {
                return mat.template topRightCorner<Eigen::Dynamic, ColSizeType::Value>(row_size, ColSizeType::Value);
            }
        }
        else if constexpr (!std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.topRightCorner(row_size.value(), col_size);
            }
            else
            {
                return mat.template topRightCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size);
            }
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic && ColSizeType::IsDynamic)
            {
                return mat.topRightCorner(row_size.value(), col_size.value());
            }
            else if constexpr (RowSizeType::IsDynamic && ColSizeType::IsFixed)
            {
                return mat.template topRightCorner<Eigen::Dynamic, ColSizeType::Value>(row_size.value(), ColSizeType::Value);
            }
            else if constexpr (RowSizeType::IsFixed && ColSizeType::IsDynamic)
            {
                return mat.template topRightCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size.value());
            }
            else
            {
                return mat.template topRightCorner<RowSizeType::Value, ColSizeType::Value>();
            }
        }
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto topRightCorner(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template topRightCorner<compile_time_rows, compile_time_cols>();
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto topRightCorner(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template topRightCorner<compile_time_rows, compile_time_cols>();
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto topRightCorner(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return topRightCornerImpl(mat, row_size, col_size);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto topRightCorner(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return topRightCornerImpl(mat, row_size, col_size);
    }

    // Bottom-left corner accessors
    template <typename MatType, typename RowSizeType, typename ColSizeType>
    static auto bottomLeftCornerImpl(MatType &&mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            return mat.bottomLeftCorner(row_size, col_size);
        }
        else if constexpr (std::is_integral_v<RowSizeType> && !std::is_integral_v<ColSizeType>)
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.bottomLeftCorner(row_size, col_size.value());
            }
            else
            {
                return mat.template bottomLeftCorner<Eigen::Dynamic, ColSizeType::Value>(row_size, ColSizeType::Value);
            }
        }
        else if constexpr (!std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.bottomLeftCorner(row_size.value(), col_size);
            }
            else
            {
                return mat.template bottomLeftCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size);
            }
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic && ColSizeType::IsDynamic)
            {
                return mat.bottomLeftCorner(row_size.value(), col_size.value());
            }
            else if constexpr (RowSizeType::IsDynamic && ColSizeType::IsFixed)
            {
                return mat.template bottomLeftCorner<Eigen::Dynamic, ColSizeType::Value>(row_size.value(), ColSizeType::Value);
            }
            else if constexpr (RowSizeType::IsFixed && ColSizeType::IsDynamic)
            {
                return mat.template bottomLeftCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size.value());
            }
            else
            {
                return mat.template bottomLeftCorner<RowSizeType::Value, ColSizeType::Value>();
            }
        }
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto bottomLeftCorner(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template bottomLeftCorner<compile_time_rows, compile_time_cols>();
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto bottomLeftCorner(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template bottomLeftCorner<compile_time_rows, compile_time_cols>();
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto bottomLeftCorner(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return bottomLeftCornerImpl(mat, row_size, col_size);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto bottomLeftCorner(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return bottomLeftCornerImpl(mat, row_size, col_size);
    }

    // Bottom-right corner accessors
    template <typename MatType, typename RowSizeType, typename ColSizeType>
    static auto bottomRightCornerImpl(MatType &&mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            return mat.bottomRightCorner(row_size, col_size);
        }
        else if constexpr (std::is_integral_v<RowSizeType> && !std::is_integral_v<ColSizeType>)
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.bottomRightCorner(row_size, col_size.value());
            }
            else
            {
                return mat.template bottomRightCorner<Eigen::Dynamic, ColSizeType::Value>(row_size, ColSizeType::Value);
            }
        }
        else if constexpr (!std::is_integral_v<RowSizeType> && std::is_integral_v<ColSizeType>)
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.bottomRightCorner(row_size.value(), col_size);
            }
            else
            {
                return mat.template bottomRightCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size);
            }
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic && ColSizeType::IsDynamic)
            {
                return mat.bottomRightCorner(row_size.value(), col_size.value());
            }
            else if constexpr (RowSizeType::IsDynamic && ColSizeType::IsFixed)
            {
                return mat.template bottomRightCorner<Eigen::Dynamic, ColSizeType::Value>(row_size.value(), ColSizeType::Value);
            }
            else if constexpr (RowSizeType::IsFixed && ColSizeType::IsDynamic)
            {
                return mat.template bottomRightCorner<RowSizeType::Value, Eigen::Dynamic>(RowSizeType::Value, col_size.value());
            }
            else
            {
                return mat.template bottomRightCorner<RowSizeType::Value, ColSizeType::Value>();
            }
        }
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto bottomRightCorner(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template bottomRightCorner<compile_time_rows, compile_time_cols>();
    }

    template <auto RowSize, auto ColSize, typename Derived>
    static auto bottomRightCorner(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template bottomRightCorner<compile_time_rows, compile_time_cols>();
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto bottomRightCorner(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return bottomRightCornerImpl(mat, row_size, col_size);
    }

    template <typename Derived, typename RowSizeType, typename ColSizeType>
    static auto bottomRightCorner(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size, const ColSizeType &col_size)
    {
        return bottomRightCornerImpl(mat, row_size, col_size);
    }

    // Row and column accessors
    template <typename MatType, typename RowSizeType>
    static auto topRowsImpl(MatType &&mat, const RowSizeType &row_size)
    {
        if constexpr (std::is_integral_v<RowSizeType>)
        {
            return mat.topRows(row_size);
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.topRows(row_size.value());
            }
            else
            {
                return mat.template topRows<RowSizeType::Value>();
            }
        }
    }

    template <auto RowSize, typename Derived>
    static auto topRows(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        return mat.template topRows<compile_time_rows>();
    }

    template <auto RowSize, typename Derived>
    static auto topRows(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        return mat.template topRows<compile_time_rows>();
    }

    template <typename Derived, typename RowSizeType>
    static auto topRows(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size)
    {
        return topRowsImpl(mat, row_size);
    }

    template <typename Derived, typename RowSizeType>
    static auto topRows(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size)
    {
        return topRowsImpl(mat, row_size);
    }

    template <typename MatType, typename RowSizeType>
    static auto bottomRowsImpl(MatType &&mat, const RowSizeType &row_size)
    {
        if constexpr (std::is_integral_v<RowSizeType>)
        {
            return mat.bottomRows(row_size);
        }
        else
        {
            if constexpr (RowSizeType::IsDynamic)
            {
                return mat.bottomRows(row_size.value());
            }
            else
            {
                return mat.template bottomRows<RowSizeType::Value>();
            }
        }
    }

    template <auto RowSize, typename Derived>
    static auto bottomRows(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        return mat.template bottomRows<compile_time_rows>();
    }

    template <auto RowSize, typename Derived>
    static auto bottomRows(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_rows = extract_compile_time_value<RowSize>::Value;
        return mat.template bottomRows<compile_time_rows>();
    }

    template <typename Derived, typename RowSizeType>
    static auto bottomRows(Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size)
    {
        return bottomRowsImpl(mat, row_size);
    }

    template <typename Derived, typename RowSizeType>
    static auto bottomRows(const Eigen::MatrixBase<Derived> &mat, const RowSizeType &row_size)
    {
        return bottomRowsImpl(mat, row_size);
    }

    template <typename MatType, typename ColSizeType>
    static auto leftColsImpl(MatType &&mat, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<ColSizeType>)
        {
            return mat.leftCols(col_size);
        }
        else
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.leftCols(col_size.value());
            }
            else
            {
                return mat.template leftCols<ColSizeType::Value>();
            }
        }
    }

    template <auto ColSize, typename Derived>
    static auto leftCols(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template leftCols<compile_time_cols>();
    }

    template <auto ColSize, typename Derived>
    static auto leftCols(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template leftCols<compile_time_cols>();
    }

    template <typename Derived, typename ColSizeType>
    static auto leftCols(Eigen::MatrixBase<Derived> &mat, const ColSizeType &col_size)
    {
        return leftColsImpl(mat, col_size);
    }

    template <typename Derived, typename ColSizeType>
    static auto leftCols(const Eigen::MatrixBase<Derived> &mat, const ColSizeType &col_size)
    {
        return leftColsImpl(mat, col_size);
    }

    template <typename MatType, typename ColSizeType>
    static auto rightColsImpl(MatType &&mat, const ColSizeType &col_size)
    {
        if constexpr (std::is_integral_v<ColSizeType>)
        {
            return mat.rightCols(col_size);
        }
        else
        {
            if constexpr (ColSizeType::IsDynamic)
            {
                return mat.rightCols(col_size.value());
            }
            else
            {
                return mat.template rightCols<ColSizeType::Value>();
            }
        }
    }

    template <auto ColSize, typename Derived>
    static auto rightCols(Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template rightCols<compile_time_cols>();
    }

    template <auto ColSize, typename Derived>
    static auto rightCols(const Eigen::MatrixBase<Derived> &mat)
    {
        constexpr int compile_time_cols = extract_compile_time_value<ColSize>::Value;
        return mat.template rightCols<compile_time_cols>();
    }

    template <typename Derived, typename ColSizeType>
    static auto rightCols(Eigen::MatrixBase<Derived> &mat, const ColSizeType &col_size)
    {
        return rightColsImpl(mat, col_size);
    }

    template <typename Derived, typename ColSizeType>
    static auto rightCols(const Eigen::MatrixBase<Derived> &mat, const ColSizeType &col_size)
    {
        return rightColsImpl(mat, col_size);
    }

} // namespace galileo

#endif // __galileo_common_meta_eigen_hpp__
