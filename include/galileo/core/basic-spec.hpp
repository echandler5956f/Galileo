#ifndef __galileo_core_basic_spec_hpp__
#define __galileo_core_basic_spec_hpp__

#include "galileo/core/fwd.hpp"

// Macros to import the types and constants from a basic spec
#define GALILEO_BASIC_SPEC_META_TYPEDEF(BasicSpec)                 \
    using RobotModel_t = typename BasicSpec::RobotModel_t;         \
    using RobotData_t = typename BasicSpec::RobotData_t;           \
    using State_t = typename BasicSpec::State_t;                   \
    using ActuationMeta_t = typename BasicSpec::ActuationMeta_t;   \
    using ActuationModel_t = typename BasicSpec::ActuationModel_t; \
    using ActuationData_t = typename BasicSpec::ActuationData_t;

#define GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(BasicSpec) \
    using VarScalar = typename BasicSpec::VarScalar;  \
    using NumScalar = typename BasicSpec::NumScalar;  \
    static constexpr int Options = BasicSpec::Options;

#define GALILEO_BASIC_SPEC_CONSTANTS_TYPEDEF(BasicSpec) \
    static constexpr int NQb = BasicSpec::NQb;          \
    static constexpr int NQj = BasicSpec::NQj;          \
    static constexpr int NVb = BasicSpec::NVb;          \
    static constexpr int NVj = BasicSpec::NVj;          \
    static constexpr int NRotors = BasicSpec::NRotors;  \
    static constexpr int NQ = BasicSpec::NQ;            \
    static constexpr int NV = BasicSpec::NV;            \
    static constexpr int NX = BasicSpec::NX;            \
    static constexpr int NDX = BasicSpec::NDX;          \
    static constexpr int NUa = BasicSpec::NUa;

#define GALILEO_BASIC_SPEC_EIGEN_TYPES_TYPEDEF(BasicSpec)    \
    using VectorNqb_t = typename BasicSpec::VectorNqb_t;     \
    using VectorNqj_t = typename BasicSpec::VectorNqj_t;     \
    using VectorNvb_t = typename BasicSpec::VectorNvb_t;     \
    using VectorNvj_t = typename BasicSpec::VectorNvj_t;     \
    using VectorNx_t = typename BasicSpec::VectorNx_t;       \
    using VectorNua_t = typename BasicSpec::VectorNua_t;     \
    using VectorNdx_t = typename BasicSpec::VectorNdx_t;     \
    using VectorNq_t = typename BasicSpec::VectorNq_t;       \
    using VectorNv_t = typename BasicSpec::VectorNv_t;       \
    using MatrixNx_t = typename BasicSpec::MatrixNx_t;       \
    using MatrixNua_t = typename BasicSpec::MatrixNua_t;     \
    using MatrixNdx_t = typename BasicSpec::MatrixNdx_t;     \
    using MatrixNq_t = typename BasicSpec::MatrixNq_t;       \
    using MatrixNv_t = typename BasicSpec::MatrixNv_t;       \
    using MatrixNvNdx_t = typename BasicSpec::MatrixNvNdx_t; \
    using MatrixNvNua_t = typename BasicSpec::MatrixNvNua_t; \
    using MatrixNuaNv_t = typename BasicSpec::MatrixNuaNv_t; \
    using MatrixNdxNua_t = typename BasicSpec::MatrixNdxNua_t;

#define GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BasicSpec) \
    GALILEO_BASIC_SPEC_META_TYPEDEF(BasicSpec);      \
    GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(BasicSpec);   \
    GALILEO_BASIC_SPEC_CONSTANTS_TYPEDEF(BasicSpec); \
    GALILEO_BASIC_SPEC_EIGEN_TYPES_TYPEDEF(BasicSpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Defines the basic types and constants used in the core library. */
    /* ---------------------------------------------------------------- */
    template <typename _VarScalar,
              typename _NumScalar,
              int _Options,
              int _NQb,
              int _NQj,
              int _NVb,
              int _NVj,
              int _NRotors,
              template <typename> class StateTpl,
              template <typename> class ActuationTpl>
    struct BasicSpecTpl
    {
        using BasicSpec = BasicSpecTpl<_VarScalar, _NumScalar, _Options, _NQb, _NQj, _NVb, _NVj, _NRotors, StateTpl, ActuationTpl>;

        /* ---------------------------------------------------------------- */
        /* Scalar types and Eigen Matrix storage order */
        /* ---------------------------------------------------------------- */
        using VarScalar = _VarScalar;            // Scalar type for variables (for AD)
        using NumScalar = _NumScalar;            // Scalar type for numerics (i.e., bounds, times, etc.)
        static constexpr int Options = _Options; // Eigen storage order

        /* ---------------------------------------------------------------- */
        /* Compile-time constants */
        /* ---------------------------------------------------------------- */
        static constexpr int NQb = _NQb;         // Dimension of floating base generalized coordinates
        static constexpr int NQj = _NQj;         // Dimension of joint generalized coordinates
        static constexpr int NVb = _NVb;         // Dimension of floating base generalized velocities
        static constexpr int NVj = _NVj;         // Dimension of joint generalized velocities
        static constexpr int NRotors = _NRotors; // Number of rotors attached to the floating base

        static constexpr int NQ = NQb + NQj; // Dimension of generalized coordinates
        static constexpr int NV = NVb + NVj; // Dimension of generalized velocities

        static constexpr int NX = NQ + NV;  // State dimension
        static constexpr int NDX = NV + NV; // State tangent space dimension

        static constexpr int NUa = NV - NVb + NRotors; // Dimension of actuated torque inputs

        /* ---------------------------------------------------------------- */
        /* Fixed-size Eigen types */
        /* ---------------------------------------------------------------- */
        using VectorNqb_t = Eigen::Matrix<VarScalar, NQb, 1, Options>;
        using VectorNqj_t = Eigen::Matrix<VarScalar, NQj, 1, Options>;
        using VectorNvb_t = Eigen::Matrix<VarScalar, NVb, 1, Options>;
        using VectorNvj_t = Eigen::Matrix<VarScalar, NVj, 1, Options>;

        using VectorNx_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
        using VectorNua_t = Eigen::Matrix<VarScalar, NUa, 1, Options>;
        using VectorNdx_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
        using VectorNq_t = Eigen::Matrix<VarScalar, NQ, 1, Options>;
        using VectorNv_t = Eigen::Matrix<VarScalar, NV, 1, Options>;

        using MatrixNx_t = Eigen::Matrix<VarScalar, NX, NX, Options>;
        using MatrixNua_t = Eigen::Matrix<VarScalar, NUa, NUa, Options>;
        using MatrixNdx_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
        using MatrixNq_t = Eigen::Matrix<VarScalar, NQ, NQ, Options>;
        using MatrixNv_t = Eigen::Matrix<VarScalar, NV, NV, Options>;

        using MatrixNvNdx_t = Eigen::Matrix<VarScalar, NV, NDX, Options>;
        using MatrixNvNua_t = Eigen::Matrix<VarScalar, NV, NUa, Options>;
        using MatrixNuaNv_t = Eigen::Matrix<VarScalar, NUa, NV, Options>;
        using MatrixNdxNua_t = Eigen::Matrix<VarScalar, NDX, NUa, Options>;

        /* ---------------------------------------------------------------- */
        /* Template types */
        /* ---------------------------------------------------------------- */
        using RobotModel_t = pinocchio::ModelTpl<VarScalar, Options>;
        using RobotData_t = pinocchio::DataTpl<VarScalar, Options>;

        using State_t = StateTpl<BasicSpec>;

        using ActuationMeta_t = ActuationTpl<BasicSpec>;
        using ActuationModel_t = typename traits<ActuationMeta_t>::Model_t;
        using ActuationData_t = typename traits<ActuationMeta_t>::Data_t;
    };

} // namespace galileo

#endif // __galileo_core_basic_spec_hpp__
