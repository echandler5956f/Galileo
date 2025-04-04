#ifndef __galileo_core_basic_spec_hpp__
#define __galileo_core_basic_spec_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    namespace core
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
                  template <typename> class RobotTpl,
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
            using Vector2_t = Eigen::Matrix<VarScalar, 2, 1, Options>;
            using Vector3_t = Eigen::Matrix<VarScalar, 3, 1, Options>;
            using Vector4_t = Eigen::Matrix<VarScalar, 4, 1, Options>;
            using Vector6_t = Eigen::Matrix<VarScalar, 6, 1, Options>;

            using VectorNqb_t = Eigen::Matrix<VarScalar, NQb, 1, Options>;
            using VectorNqj_t = Eigen::Matrix<VarScalar, NQj, 1, Options>;
            using VectorNvb_t = Eigen::Matrix<VarScalar, NVb, 1, Options>;
            using VectorNvj_t = Eigen::Matrix<VarScalar, NVj, 1, Options>;

            using VectorNx_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
            using VectorNua_t = Eigen::Matrix<VarScalar, NUa, 1, Options>;
            using VectorNdx_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using VectorNq_t = Eigen::Matrix<VarScalar, NQ, 1, Options>;
            using VectorNv_t = Eigen::Matrix<VarScalar, NV, 1, Options>;

            using Matrix2_t = Eigen::Matrix<VarScalar, 2, 2, Options>;
            using Matrix3_t = Eigen::Matrix<VarScalar, 3, 3, Options>;
            using Matrix4_t = Eigen::Matrix<VarScalar, 4, 4, Options>;
            using Matrix6_t = Eigen::Matrix<VarScalar, 6, 6, Options>;

            using Matrix3Ndx_t = Eigen::Matrix<VarScalar, 3, NDX, Options>;
            using Matrix3Nv_t = Eigen::Matrix<VarScalar, 3, NV, Options>;
            using Matrix6Ndx_t = Eigen::Matrix<VarScalar, 6, NDX, Options>;
            using Matrix6Nv_t = Eigen::Matrix<VarScalar, 6, NV, Options>;

            using MatrixNx_t = Eigen::Matrix<VarScalar, NX, NX, Options>;
            using MatrixNua_t = Eigen::Matrix<VarScalar, NUa, NUa, Options>;
            using MatrixNdx_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
            using MatrixNq_t = Eigen::Matrix<VarScalar, NQ, NQ, Options>;
            using MatrixNv_t = Eigen::Matrix<VarScalar, NV, NV, Options>;

            using MatrixNdxNua_t = Eigen::Matrix<VarScalar, NDX, NUa, Options>;
            using MatrixNuaNv_t = Eigen::Matrix<VarScalar, NUa, NV, Options>;
            using MatrixNvNdx_t = Eigen::Matrix<VarScalar, NV, NDX, Options>;
            using MatrixNvNua_t = Eigen::Matrix<VarScalar, NV, NUa, Options>;

            /* ---------------------------------------------------------------- */
            /* Template types */
            /* ---------------------------------------------------------------- */

            using RobotMeta_t = RobotTpl<BasicSpec>;
            using RobotModel_t = typename RobotMeta_t::Model;
            using RobotData_t = typename RobotMeta_t::Data;

            using State_t = StateTpl<BasicSpec>;

            using ActuationMeta_t = ActuationTpl<BasicSpec>;
            using ActuationModel_t = typename ActuationMeta_t::Model;
            using ActuationData_t = typename ActuationMeta_t::Data;
        };

    } // namespace core

} // namespace galileo

#endif // __galileo_core_basic_spec_hpp__
