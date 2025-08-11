#ifndef __galileo_multibody_robot_spec_hpp__
#define __galileo_multibody_robot_spec_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/force.hpp>
#include <pinocchio/spatial/motion.hpp>
#include <pinocchio/spatial/se3.hpp>

#include "galileo/multibody/robot-holder.hpp"

#include "galileo/core/basic-spec.hpp"

// Macros to import the types and constants from a robot spec
#define GALILEO_ROBOT_SPEC_META_TYPEDEF(RobotSpec)                 \
    using State_t = typename RobotSpec::State_t;                   \
    using ActuationMeta_t = typename RobotSpec::ActuationMeta_t;   \
    using ActuationModel_t = typename RobotSpec::ActuationModel_t; \
    using ActuationData_t = typename RobotSpec::ActuationData_t;

#define GALILEO_ROBOT_SPEC_PINOCCHIO_TYPES_TYPEDEF(RobotSpec)      \
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
    using Matrix6Nv_t = typename RobotSpec::Matrix6Nv_t;               \
    using Matrix6Ndx_t = typename RobotSpec::Matrix6Ndx_t;

#define GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(RobotSpec)       \
    GALILEO_ROBOT_SPEC_META_TYPEDEF(RobotSpec);            \
    GALILEO_ROBOT_SPEC_PINOCCHIO_TYPES_TYPEDEF(RobotSpec); \
    GALILEO_ROBOT_SPEC_SCALARS_TYPEDEF(RobotSpec);         \
    GALILEO_ROBOT_SPEC_CONSTANTS_TYPEDEF(RobotSpec);       \
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
        : public RobotHolderTpl<_NQb, _NQj, _NVb, _NVj, _NRotors>
    {
        using RS = RobotSpecTpl<BasicSpec, _NQb, _NQj, _NVb, _NVj, _NRotors, StateTpl, ActuationTpl>;
        using Base = RobotHolderTpl<_NQb, _NQj, _NVb, _NVj, _NRotors>;

        /* ---------------------------------------------------------------- */
        /* Forward the dimension types from the robot holder */
        /* ---------------------------------------------------------------- */
        using Base::DimNDX_t;
        using Base::DimNQ_t;
        using Base::DimNQb_t;
        using Base::DimNQj_t;
        using Base::DimNRotors_t;
        using Base::DimNUa_t;
        using Base::DimNV_t;
        using Base::DimNVb_t;
        using Base::DimNVj_t;
        using Base::DimNX_t;

        /* ------------------------------------------------------------------- */
        /* Forward the associated compile-time constants from the robot holder */
        /* ------------------------------------------------------------------- */
        static constexpr int NQb = Base::NQb;
        static constexpr int NQj = Base::NQj;
        static constexpr int NVb = Base::NVb;
        static constexpr int NVj = Base::NVj;
        static constexpr int NRotors = Base::NRotors;
        static constexpr int NQ = Base::NQ;
        static constexpr int NV = Base::NV;
        static constexpr int NX = Base::NX;
        static constexpr int NDX = Base::NDX;
        static constexpr int NUa = Base::NUa;

        /* ---------------------------------------------------------------- */
        /* Import the basic spec types and constants */
        /* ---------------------------------------------------------------- */
        using BS = BasicSpec;
        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

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

        using Base::get_nqb;
        using Base::get_nqb_dim;

        using Base::get_nqj;
        using Base::get_nqj_dim;

        using Base::get_nvb;
        using Base::get_nvb_dim;

        using Base::get_nvj;
        using Base::get_nvj_dim;

        using Base::get_nrotors;
        using Base::get_nrotors_dim;

        using Base::get_nq;
        using Base::get_nq_dim;

        using Base::get_nv;
        using Base::get_nv_dim;

        using Base::get_nx;
        using Base::get_nx_dim;

        using Base::get_ndx;
        using Base::get_ndx_dim;

        using Base::get_nua;
        using Base::get_nua_dim;

        // Constructor to properly initialize compound dimensions through the base class
        RobotSpecTpl()
            : Base()
        {
        }

        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            os << indent << "BasicSpec: {";
            BS{}.display(os, "");
            os << "}\n";

            os << indent << "RobotHolder: {\n";
            static_cast<const Base &>(*this).display(os, indent + "  ");
            os << indent << "}\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const RobotSpecTpl &rs)
        {
            os << "RobotSpec: {\n";
            rs.display(os, "  ");
            os << "}";
            return os;
        }

    };

} // namespace galileo

#endif // __galileo_multibody_robot_spec_hpp__
