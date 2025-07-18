#ifndef __galileo_multibody_robot_spec_hpp__
#define __galileo_multibody_robot_spec_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/force.hpp>
#include <pinocchio/spatial/motion.hpp>
#include <pinocchio/spatial/se3.hpp>

#include "galileo/core/fwd.hpp"

#include "galileo/core/basic-spec.hpp"

// Macros to import the types and constants from a robot spec
#define GALILEO_ROBOT_SPEC_META_TYPEDEF(RobotSpec)                 \
    using State_t = typename RobotSpec::State_t;                   \
    using ActuationMeta_t = typename RobotSpec::ActuationMeta_t;   \
    using ActuationModel_t = typename RobotSpec::ActuationModel_t; \
    using ActuationData_t = typename RobotSpec::ActuationData_t;

#define GALILEO_ROBOT_SPEC_PINOCCIO_TYPES_TYPEDEF(RobotSpec)       \
    using RobotModel_t = typename RobotSpec::RobotModel_t;         \
    using RobotData_t = typename RobotSpec::RobotData_t;           \
    using FrameIndex_t = typename RobotSpec::FrameIndex_t;         \
    using ReferenceFrame_t = typename RobotSpec::ReferenceFrame_t; \
    using SE3_t = typename RobotSpec::SE3_t;                       \
    using Motion_t = typename RobotSpec::Motion_t;                 \
    using Force_t = typename RobotSpec::Force_t;                   \
    using ActionMatrix_t = typename RobotSpec::ActionMatrix_t;

#define GALILEO_ROBOT_SPEC_SCALARS_TYPEDEF(RobotSpec) \
    GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(RobotSpec::BS);

#define GALILEO_ROBOT_SPEC_CONSTANTS_TYPEDEF(RobotSpec)    \
    static constexpr int NQb = RobotSpec::NQb;             \
    static constexpr int NQj = RobotSpec::NQj;             \
    static constexpr int NVb = RobotSpec::NVb;             \
    static constexpr int NVj = RobotSpec::NVj;             \
    static constexpr int NRotors = RobotSpec::NRotors;     \
    static constexpr int NQ = RobotSpec::NQ;               \
    static constexpr int NV = RobotSpec::NV;               \
    static constexpr int NX = RobotSpec::NX;               \
    static constexpr int NDX = RobotSpec::NDX;             \
    static constexpr int NUa = RobotSpec::NUa;             \
    using DimNQb_t = typename RobotSpec::DimNQb_t;         \
    using DimNQj_t = typename RobotSpec::DimNQj_t;         \
    using DimNVb_t = typename RobotSpec::DimNVb_t;         \
    using DimNVj_t = typename RobotSpec::DimNVj_t;         \
    using DimNRotors_t = typename RobotSpec::DimNRotors_t; \
    using DimNQ_t = typename RobotSpec::DimNQ_t;           \
    using DimNV_t = typename RobotSpec::DimNV_t;           \
    using DimNX_t = typename RobotSpec::DimNX_t;           \
    using DimNDX_t = typename RobotSpec::DimNDX_t;         \
    using DimNUa_t = typename RobotSpec::DimNUa_t;

#define GALILEO_ROBOT_SPEC_EIGEN_TYPES_TYPEDEF(RobotSpec)              \
    GALILEO_BASIC_SPEC_FIXED_SIZE_EIGEN_TYPES_TYPEDEF(RobotSpec::BS)   \
    GALILEO_BASIC_SPEC_DYNAMIC_SIZE_EIGEN_TYPES_TYPEDEF(RobotSpec::BS) \
    using VectorNqb_t = typename RobotSpec::VectorNqb_t;               \
    using VectorNqj_t = typename RobotSpec::VectorNqj_t;               \
    using VectorNvb_t = typename RobotSpec::VectorNvb_t;               \
    using VectorNvj_t = typename RobotSpec::VectorNvj_t;               \
    using VectorNx_t = typename RobotSpec::VectorNx_t;                 \
    using VectorNua_t = typename RobotSpec::VectorNua_t;               \
    using VectorNdx_t = typename RobotSpec::VectorNdx_t;               \
    using VectorNq_t = typename RobotSpec::VectorNq_t;                 \
    using VectorNv_t = typename RobotSpec::VectorNv_t;                 \
    using MatrixNx_t = typename RobotSpec::MatrixNx_t;                 \
    using MatrixNua_t = typename RobotSpec::MatrixNua_t;               \
    using MatrixNdx_t = typename RobotSpec::MatrixNdx_t;               \
    using MatrixNq_t = typename RobotSpec::MatrixNq_t;                 \
    using MatrixNv_t = typename RobotSpec::MatrixNv_t;                 \
    using MatrixNvNdx_t = typename RobotSpec::MatrixNvNdx_t;           \
    using MatrixNvNua_t = typename RobotSpec::MatrixNvNua_t;           \
    using MatrixNuaNv_t = typename RobotSpec::MatrixNuaNv_t;           \
    using MatrixNdxNua_t = typename RobotSpec::MatrixNdxNua_t;         \
    using MatrixNuaNdx_t = typename RobotSpec::MatrixNuaNdx_t;         \
    using MatrixNv6_t = typename RobotSpec::MatrixNv6_t;               \
    using Matrix6Nv_t = typename RobotSpec::Matrix6Nv_t;

#define GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RobotSpec)      \
    GALILEO_ROBOT_SPEC_META_TYPEDEF(RobotSpec);           \
    GALILEO_ROBOT_SPEC_PINOCCIO_TYPES_TYPEDEF(RobotSpec); \
    GALILEO_ROBOT_SPEC_SCALARS_TYPEDEF(RobotSpec);        \
    GALILEO_ROBOT_SPEC_CONSTANTS_TYPEDEF(RobotSpec);      \
    GALILEO_ROBOT_SPEC_EIGEN_TYPES_TYPEDEF(RobotSpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Defines the robot-specific types and constants used in the core library. */
    /* ---------------------------------------------------------------- */
    template <typename BasicSpec,
              int _NQb,
              int _NQj,
              int _NVb,
              int _NVj,
              int _NRotors,
              template <typename> class StateTpl,
              template <typename> class ActuationTpl>
    struct RobotSpecTpl
    {
        using BS = BasicSpec;
        using RS = RobotSpecTpl<BS, _NQb, _NQj, _NVb, _NVj, _NRotors, StateTpl, ActuationTpl>;

        // Import the basic spec types and constants
        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

        /* ---------------------------------------------------------------- */
        /* Dimension types */
        /* ---------------------------------------------------------------- */
        using DimNQb_t = DimensionTpl<_NQb>;
        using DimNQj_t = DimensionTpl<_NQj>;
        using DimNVb_t = DimensionTpl<_NVb>;
        using DimNVj_t = DimensionTpl<_NVj>;
        using DimNRotors_t = DimensionTpl<_NRotors>;
        using DimNQ_t = decltype(DimNQb_t{} + DimNQj_t{});
        using DimNV_t = decltype(DimNVb_t{} + DimNVj_t{});
        using DimNX_t = decltype(DimNQ_t{} + DimNV_t{});
        using DimNDX_t = decltype(DimNV_t{} + DimNV_t{});
        using DimNUa_t = decltype(DimNVj_t{} + DimNRotors_t{});

        /* ---------------------------------------------------------------- */
        /* Compile-time constants */
        /* ---------------------------------------------------------------- */
        static constexpr int NQb = DimNQb_t::Value;         // Dimension of floating base generalized coordinates
        static constexpr int NQj = DimNQj_t::Value;         // Dimension of joint generalized coordinates
        static constexpr int NVb = DimNVb_t::Value;         // Dimension of floating base generalized velocities
        static constexpr int NVj = DimNVj_t::Value;         // Dimension of joint generalized velocities
        static constexpr int NRotors = DimNRotors_t::Value; // Number of rotors attached to the floating base

        static constexpr int NQ = DimNQ_t::Value; // Dimension of generalized coordinates
        static constexpr int NV = DimNV_t::Value; // Dimension of generalized velocities

        static constexpr int NX = DimNX_t::Value;   // State dimension
        static constexpr int NDX = DimNDX_t::Value; // State tangent space dimension

        static constexpr int NUa = DimNUa_t::Value; // Dimension of actuated torque inputs

        /* ---------------------------------------------------------------- */
        /* Fixed-size Eigen types */
        /* ---------------------------------------------------------------- */
        using VectorNqb_t = Eigen::GMatrix<VarScalar, NQb, 1, Options>;
        using VectorNqj_t = Eigen::GMatrix<VarScalar, NQj, 1, Options>;
        using VectorNvb_t = Eigen::GMatrix<VarScalar, NVb, 1, Options>;
        using VectorNvj_t = Eigen::GMatrix<VarScalar, NVj, 1, Options>;

        using VectorNx_t = Eigen::GMatrix<VarScalar, NX, 1, Options>;
        using VectorNua_t = Eigen::GMatrix<VarScalar, NUa, 1, Options>;
        using VectorNdx_t = Eigen::GMatrix<VarScalar, NDX, 1, Options>;
        using VectorNq_t = Eigen::GMatrix<VarScalar, NQ, 1, Options>;
        using VectorNv_t = Eigen::GMatrix<VarScalar, NV, 1, Options>;

        using MatrixNx_t = Eigen::GMatrix<VarScalar, NX, NX, Options>;
        using MatrixNua_t = Eigen::GMatrix<VarScalar, NUa, NUa, Options>;
        using MatrixNdx_t = Eigen::GMatrix<VarScalar, NDX, NDX, Options>;
        using MatrixNq_t = Eigen::GMatrix<VarScalar, NQ, NQ, Options>;
        using MatrixNv_t = Eigen::GMatrix<VarScalar, NV, NV, Options>;

        using MatrixNvNdx_t = Eigen::GMatrix<VarScalar, NV, NDX, Options>;
        using MatrixNvNua_t = Eigen::GMatrix<VarScalar, NV, NUa, Options>;
        using MatrixNuaNv_t = Eigen::GMatrix<VarScalar, NUa, NV, Options>;
        using MatrixNdxNua_t = Eigen::GMatrix<VarScalar, NDX, NUa, Options>;
        using MatrixNuaNdx_t = Eigen::GMatrix<VarScalar, NUa, NDX, Options>;

        using MatrixNv6_t = Eigen::GMatrix<VarScalar, NV, 6, Options>;
        using Matrix6Nv_t = Eigen::GMatrix<VarScalar, 6, NV, Options>;
        using Matrix6Ndx_t = Eigen::GMatrix<VarScalar, 6, NDX, Options>;

        /* ---------------------------------------------------------------- */
        /* Template types */
        /* ---------------------------------------------------------------- */
        using State_t = StateTpl<RS>;

        using ActuationMeta_t = ActuationTpl<RS>;
        using ActuationModel_t = typename traits<ActuationMeta_t>::Model_t;
        using ActuationData_t = typename traits<ActuationMeta_t>::Data_t;

        /* ---------------------------------------------------------------- */
        /* Pinocchio types */
        /* ---------------------------------------------------------------- */
        using RobotModel_t = pinocchio::ModelTpl<VarScalar, Options>;
        using RobotData_t = pinocchio::DataTpl<VarScalar, Options>;
        using FrameIndex_t = pinocchio::FrameIndex;
        using ReferenceFrame_t = pinocchio::ReferenceFrame;
        using SE3_t = pinocchio::SE3Tpl<VarScalar, Options>;
        using Motion_t = pinocchio::MotionTpl<VarScalar, Options>;
        using Force_t = pinocchio::ForceTpl<VarScalar, Options>;
        using ActionMatrix_t = typename SE3_t::ActionMatrixType;

        /* ---------------------------------------------------------------- */
        /* Accessors for the RobotSpec dimensions */
        /* ---------------------------------------------------------------- */

        const DimNQb_t &get_nqb_dim() const
        {
            return nqb_dim_;
        }

        const int get_nqb() const
        {
            return nqb_dim_.value();
        }

        const DimNQj_t &get_nqj_dim() const
        {
            return nqj_dim_;
        }

        const int get_nqj() const
        {
            return nqj_dim_.value();
        }

        const DimNVb_t &get_nvb_dim() const
        {
            return nvb_dim_;
        }

        const int get_nvb() const
        {
            return nvb_dim_.value();
        }

        const DimNVj_t &get_nvj_dim() const
        {
            return nvj_dim_;
        }

        const int get_nvj() const
        {
            return nvj_dim_.value();
        }

        const DimNRotors_t &get_nrotors_dim() const
        {
            return nrotors_dim_;
        }

        const int get_nrotors() const
        {
            return nrotors_dim_.value();
        }

        const DimNQ_t &get_nq_dim() const
        {
            return nq_dim_;
        }

        const int get_nq() const
        {
            return nq_dim_.value();
        }

        const DimNV_t &get_nv_dim() const
        {
            return nv_dim_;
        }

        const int get_nv() const
        {
            return nv_dim_.value();
        }

        const DimNX_t &get_nx_dim() const
        {
            return nx_dim_;
        }

        const int get_nx() const
        {
            return nx_dim_.value();
        }

        const DimNDX_t &get_ndx_dim() const
        {
            return ndx_dim_;
        }

        const int get_ndx() const
        {
            return ndx_dim_.value();
        }

        const DimNUa_t &get_nua_dim() const
        {
            return nua_dim_;
        }

        const int get_nua() const
        {
            return nua_dim_.value();
        }

        // Constructor to properly initialize compound dimensions
        RobotSpecTpl()
            : nqb_dim_{}, nqj_dim_{}, nvb_dim_{}, nvj_dim_{}, nrotors_dim_{},
              nq_dim_(nqb_dim_ + nqj_dim_),
              nv_dim_(nvb_dim_ + nvj_dim_),
              nx_dim_(nq_dim_ + nv_dim_),
              ndx_dim_(nv_dim_ + nv_dim_),
              nua_dim_(nvj_dim_ + nrotors_dim_)
        {
        }

        /* ---------------------------------------------------------------- */
        /* Actual storage of dimension types */
        /* ---------------------------------------------------------------- */
        DimNQb_t nqb_dim_;
        DimNQj_t nqj_dim_;
        DimNVb_t nvb_dim_;
        DimNVj_t nvj_dim_;
        DimNRotors_t nrotors_dim_;
        DimNQ_t nq_dim_;
        DimNV_t nv_dim_;
        DimNX_t nx_dim_;
        DimNDX_t ndx_dim_;
        DimNUa_t nua_dim_;
    };

    // Helper function to validate if a robot spec is in a valid configuration at runtime
    template <typename RobotSpec>
    bool IsValidRobotSpec(const RobotSpec &rs)
    {
        bool valid_nq = (rs.get_nq() == rs.get_nqb() + rs.get_nqj());
        bool valid_nv = (rs.get_nv() == rs.get_nvb() + rs.get_nvj());
        bool valid_nx = (rs.get_nx() == rs.get_nq() + rs.get_nv());
        bool valid_ndx = (rs.get_ndx() == rs.get_nv() + rs.get_nv());
        bool valid_nua = (rs.get_nua() == rs.get_nvj() + rs.get_nrotors());
        return valid_nq && valid_nv && valid_nx && valid_ndx && valid_nua;
    }

} // namespace galileo

#endif // __galileo_multibody_robot_spec_hpp__
