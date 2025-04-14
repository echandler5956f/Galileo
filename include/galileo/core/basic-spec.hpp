#ifndef __galileo_core_basic_spec_hpp__
#define __galileo_core_basic_spec_hpp__

#include "galileo/core/fwd.hpp"

#define GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(BasicSpec) \
    using VarScalar = typename BasicSpec::VarScalar;  \
    using NumScalar = typename BasicSpec::NumScalar;  \
    static constexpr int Options = BasicSpec::Options;

#define GALILEO_BASIC_SPEC_FIXED_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec) \
    using Vector2s = typename BasicSpec::Vector2s;                   \
    using Vector3s = typename BasicSpec::Vector3s;                   \
    using Vector4s = typename BasicSpec::Vector4s;                   \
    using Vector6s = typename BasicSpec::Vector6s;                   \
    using Matrix2s = typename BasicSpec::Matrix2s;                   \
    using Matrix3s = typename BasicSpec::Matrix3s;                   \
    using Matrix46s = typename BasicSpec::Matrix46s;                 \
    using Matrix6s = typename BasicSpec::Matrix6s;                   \
    using RowVector2s = typename BasicSpec::RowVector2s;

#define GALILEO_BASIC_SPEC_DYNAMIC_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec) \
    using MatrixX3s = typename BasicSpec::MatrixX3s;                   \
    using MatrixX6s = typename BasicSpec::MatrixX6s;                   \
    using Matrix3xs = typename BasicSpec::Matrix3xs;                   \
    using Matrix6xs = typename BasicSpec::Matrix6xs;                   \
    using VectorXs = typename BasicSpec::VectorXs;                     \
    using MatrixXs = typename BasicSpec::MatrixXs;                     \
    using MatrixXsRowMajor = typename BasicSpec::MatrixXsRowMajor;     \
    using ArrayXs = typename BasicSpec::ArrayXs;                       \
    using Quaternions = typename BasicSpec::Quaternions;               \
    using DiagonalMatrixXs = typename BasicSpec::DiagonalMatrixXs;

#define GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BasicSpec)              \
    GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(BasicSpec);                \
    GALILEO_BASIC_SPEC_FIXED_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec); \
    GALILEO_BASIC_SPEC_DYNAMIC_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Defines the basic types used in the core library. */
    /* ---------------------------------------------------------------- */
    template <typename _VarScalar,
              typename _NumScalar,
              int _Options>
    struct BasicSpecTpl
    {
        using BasicSpec = BasicSpecTpl<_VarScalar, _NumScalar, _Options>;

        /* ---------------------------------------------------------------- */
        /* Scalar types and Eigen Matrix storage order */
        /* ---------------------------------------------------------------- */
        using VarScalar = _VarScalar;            // Scalar type for variables (for AD)
        using NumScalar = _NumScalar;            // Scalar type for numerics (i.e., bounds, times, etc.)
        static constexpr int Options = _Options; // Eigen storage order

        /* ---------------------------------------------------------------- */
        /* Fixed-size Eigen types */
        /* ---------------------------------------------------------------- */
        using Vector2s = Eigen::Matrix<VarScalar, 2, 1, Options>;
        using Vector3s = Eigen::Matrix<VarScalar, 3, 1, Options>;
        using Vector4s = Eigen::Matrix<VarScalar, 4, 1, Options>;
        using Vector6s = Eigen::Matrix<VarScalar, 6, 1, Options>;
        using Matrix2s = Eigen::Matrix<VarScalar, 2, 2, Options>;
        using Matrix3s = Eigen::Matrix<VarScalar, 3, 3, Options>;
        using Matrix46s = Eigen::Matrix<VarScalar, 4, 6, Options>;
        using Matrix6s = Eigen::Matrix<VarScalar, 6, 6, Options>;
        using RowVector2s = Eigen::Matrix<VarScalar, 1, 2, Options>;

        using MatrixX3s = Eigen::Matrix<VarScalar, Eigen::Dynamic, 3, Options>;
        using MatrixX6s = Eigen::Matrix<VarScalar, Eigen::Dynamic, 6, Options>;
        using Matrix3xs = Eigen::Matrix<VarScalar, 3, Eigen::Dynamic, Options>;
        using Matrix6xs = Eigen::Matrix<VarScalar, 6, Eigen::Dynamic, Options>;

        using VectorXs = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;
        using MatrixXs = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Options>;
        using MatrixXsRowMajor = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using ArrayXs = Eigen::Array<VarScalar, Eigen::Dynamic, 1>;
        using Quaternions = Eigen::Quaternion<VarScalar>;
        using DiagonalMatrixXs = Eigen::DiagonalMatrix<VarScalar, Eigen::Dynamic>;
    };

} // namespace galileo

#endif // __galileo_core_basic_spec_hpp__
