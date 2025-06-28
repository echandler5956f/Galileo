#ifndef __galileo_core_basic_spec_hpp__
#define __galileo_core_basic_spec_hpp__

#include "galileo/core/fwd.hpp"

#define GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(BasicSpec) \
    using VarScalar = typename BasicSpec::VarScalar;  \
    using NumScalar = typename BasicSpec::NumScalar;  \
    static constexpr int Options = BasicSpec::Options;

#define GALILEO_BASIC_SPEC_FIXED_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec) \
    using Vector2_t = typename BasicSpec::Vector2_t;                 \
    using Vector3_t = typename BasicSpec::Vector3_t;                 \
    using Vector4_t = typename BasicSpec::Vector4_t;                 \
    using Vector6_t = typename BasicSpec::Vector6_t;                 \
    using Matrix2_t = typename BasicSpec::Matrix2_t;                 \
    using Matrix3_t = typename BasicSpec::Matrix3_t;                 \
    using Matrix46_t = typename BasicSpec::Matrix46_t;               \
    using Matrix6_t = typename BasicSpec::Matrix6_t;                 \
    using RowVector2_t = typename BasicSpec::RowVector2_t;

#define GALILEO_BASIC_SPEC_DYNAMIC_SIZE_EIGEN_TYPES_TYPEDEF(BasicSpec) \
    using MatrixX3_t = typename BasicSpec::MatrixX3_t;                 \
    using MatrixX6_t = typename BasicSpec::MatrixX6_t;                 \
    using Matrix3X_t = typename BasicSpec::Matrix3X_t;                 \
    using Matrix6X_t = typename BasicSpec::Matrix6X_t;                 \
    using VectorX_t = typename BasicSpec::VectorX_t;                   \
    using MatrixX_t = typename BasicSpec::MatrixX_t;                   \
    using MatrixXRowMajor_t = typename BasicSpec::MatrixXRowMajor_t;   \
    using ArrayX_t = typename BasicSpec::ArrayX_t;                     \
    using Quaternion_t = typename BasicSpec::Quaternion_t;             \
    using DiagonalMatrixX_t = typename BasicSpec::DiagonalMatrixX_t;

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
        using Vector2_t = Eigen::Matrix<VarScalar, 2, 1, Options>;
        using Vector3_t = Eigen::Matrix<VarScalar, 3, 1, Options>;
        using Vector4_t = Eigen::Matrix<VarScalar, 4, 1, Options>;
        using Vector6_t = Eigen::Matrix<VarScalar, 6, 1, Options>;
        using Matrix2_t = Eigen::Matrix<VarScalar, 2, 2, Options>;
        using Matrix3_t = Eigen::Matrix<VarScalar, 3, 3, Options>;
        using Matrix46_t = Eigen::Matrix<VarScalar, 4, 6, Options>;
        using Matrix6_t = Eigen::Matrix<VarScalar, 6, 6, Options>;
        using RowVector2_t = galileo::Matrix<VarScalar, 1, 2, Options>;

        using MatrixX3_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 3, Options>;
        using MatrixX6_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 6, Options>;
        using Matrix3X_t = Eigen::Matrix<VarScalar, 3, Eigen::Dynamic, Options>;
        using Matrix6X_t = Eigen::Matrix<VarScalar, 6, Eigen::Dynamic, Options>;

        using VectorX_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;
        using MatrixX_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Options>;
        using MatrixXRowMajor_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Options | Eigen::RowMajor>;
        using ArrayX_t = Eigen::Array<VarScalar, Eigen::Dynamic, 1>;
        using Quaternion_t = Eigen::Quaternion<VarScalar>;
        using DiagonalMatrixX_t = Eigen::DiagonalMatrix<VarScalar, Eigen::Dynamic>;
    };

} // namespace galileo

#endif // __galileo_core_basic_spec_hpp__
