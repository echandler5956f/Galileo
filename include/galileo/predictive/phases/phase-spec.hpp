#ifndef __galileo_predictive_phases_phase_spec_hpp__
#define __galileo_predictive_phases_phase_spec_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include "galileo/core/basic-spec.hpp"

#include <vector>
#include <array>

#define GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec)                                 \
    using RobotMeta_t = typename PhaseSpec::RobotMeta_t;                           \
    using RobotModel_t = typename PhaseSpec::RobotModel_t;                         \
    using RobotData_t = typename PhaseSpec::RobotData_t;                           \
    using State_t = typename PhaseSpec::State_t;                                   \
    using ActuationMeta_t = typename PhaseSpec::ActuationMeta_t;                   \
    using ActuationModel_t = typename PhaseSpec::ActuationModel_t;                 \
    using ActuationData_t = typename PhaseSpec::ActuationData_t;                   \
    using ConstraintManagerMeta_t = typename PhaseSpec::ConstraintManagerMeta_t;   \
    using ConstraintCollection_t = typename PhaseSpec::ConstraintCollection_t;     \
    using ConstraintModelManager_t = typename PhaseSpec::ConstraintModelManager_t; \
    using ConstraintDataManager_t = typename PhaseSpec::ConstraintDataManager_t;   \
    using CostManagerMeta_t = typename PhaseSpec::CostManagerMeta_t;               \
    using CostCollection_t = typename PhaseSpec::CostCollection_t;                 \
    using CostModelManager_t = typename PhaseSpec::CostModelManager_t;             \
    using CostDataManager_t = typename PhaseSpec::CostDataManager_t;               \
    using NodeMeta_t = typename PhaseSpec::NodeMeta_t;                             \
    using NodeModel_t = typename PhaseSpec::NodeModel_t;                           \
    using NodeData_t = typename PhaseSpec::NodeData_t;                             \
    using NodeDataVector_t = typename PhaseSpec::NodeDataVector_t;                 \
    using ControlParamMeta_t = typename PhaseSpec::ControlParamMeta_t;             \
    using ControlParamModel_t = typename PhaseSpec::ControlParamModel_t;           \
    using ControlParamData_t = typename PhaseSpec::ControlParamData_t;             \
    using ControlParamDataVector_t = typename PhaseSpec::ControlParamDataVector_t; \
    using SegmentMeta_t = typename PhaseSpec::SegmentMeta_t;                       \
    using SegmentModel_t = typename PhaseSpec::SegmentModel_t;                     \
    using SegmentData_t = typename PhaseSpec::SegmentData_t;                       \
    using SegmentDataVector_t = typename PhaseSpec::SegmentDataVector_t;           \
    using PhaseMeta_t = typename PhaseSpec::PhaseMeta_t;                           \
    using PhaseModel_t = typename PhaseSpec::PhaseModel_t;                         \
    using PhaseData_t = typename PhaseSpec::PhaseData_t;                           \
    using PhaseDataVector_t = typename PhaseSpec::PhaseDataVector_t;

#define GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PhaseSpec) \
    using VarScalar = typename PhaseSpec::VarScalar;  \
    using NumScalar = typename PhaseSpec::NumScalar;  \
    static constexpr int Options = PhaseSpec::Options;

#define GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec) \
    static constexpr int NX = PhaseSpec::NX;            \
    static constexpr int NDX = PhaseSpec::NDX;          \
    static constexpr int NQ = PhaseSpec::NQ;            \
    static constexpr int NV = PhaseSpec::NV;            \
    static constexpr int NU = PhaseSpec::NU;            \
    static constexpr int NOrder = PhaseSpec::NOrder;    \
    static constexpr int NW = PhaseSpec::NW;            \
    static constexpr int NStages = PhaseSpec::NStages;  \
    static constexpr SegmentType SegmentType = PhaseSpec::SegmentType;

#define GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec)   \
    using XAcc_t = typename PhaseSpec::XAcc_t;             \
    using XAccx_t = typename PhaseSpec::XAccx_t;           \
    using XAccu_t = typename PhaseSpec::XAccu_t;           \
    using L_t = typename PhaseSpec::L_t;                   \
    using Lx_t = typename PhaseSpec::Lx_t;                 \
    using Lu_t = typename PhaseSpec::Lu_t;                 \
    using Lxx_t = typename PhaseSpec::Lxx_t;               \
    using Lxu_t = typename PhaseSpec::Lxu_t;               \
    using Luu_t = typename PhaseSpec::Luu_t;               \
    using H_t = typename PhaseSpec::H_t;                   \
    using Hx_t = typename PhaseSpec::Hx_t;                 \
    using Hu_t = typename PhaseSpec::Hu_t;                 \
    using H_Equality_t = typename PhaseSpec::H_Equality_t; \
    using G_t = typename PhaseSpec::G_t;                   \
    using Gx_t = typename PhaseSpec::Gx_t;                 \
    using Gu_t = typename PhaseSpec::Gu_t;                 \
    using G_Bound_t = typename PhaseSpec::G_Bound_t;

#define GALILEO_PHASE_SPEC_CONTROL_PARAM_TYPES_TYPEDEF(PhaseSpec) \
    using U_t = typename PhaseSpec::U_t;                          \
    using W_t = typename PhaseSpec::W_t;                          \
    using Uw_t = typename PhaseSpec::Uw_t;

#define GALILEO_PHASE_SPEC_SEGMENT_TYPES_TYPEDEF(PhaseSpec)              \
    using Timings_t = typename PhaseSpec::Timings_t;                     \
    using Quadrature_t = typename PhaseSpec::Quadrature_t;               \
    using StageCoefficients_t = typename PhaseSpec::StageCoefficients_t; \
    using XNext_t = typename PhaseSpec::XNext_t;                         \
    using Fx_t = typename PhaseSpec::Fx_t;                               \
    using Fw_t = typename PhaseSpec::Fw_t;                               \
    using Lw_t = typename PhaseSpec::Lw_t;                               \
    using Lxw_t = typename PhaseSpec::Lxw_t;                             \
    using Lww_t = typename PhaseSpec::Lww_t;                             \
    using Hw_t = typename PhaseSpec::Hw_t;                               \
    using Gw_t = typename PhaseSpec::Gw_t;

#define GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec)    \
    using VectorNx_t = typename PhaseSpec::VectorNx_t;       \
    using VectorNu_t = typename PhaseSpec::VectorNu_t;       \
    using VectorNdx_t = typename PhaseSpec::VectorNdx_t;     \
    using VectorNq_t = typename PhaseSpec::VectorNq_t;       \
    using VectorNv_t = typename PhaseSpec::VectorNv_t;       \
    using MatrixNx_t = typename PhaseSpec::MatrixNx_t;       \
    using MatrixNu_t = typename PhaseSpec::MatrixNu_t;       \
    using MatrixNdx_t = typename PhaseSpec::MatrixNdx_t;     \
    using MatrixNv_t = typename PhaseSpec::MatrixNv_t;       \
    using MatrixNvNw_t = typename PhaseSpec::MatrixNvNw_t;   \
    using MatrixNvNdx_t = typename PhaseSpec::MatrixNvNdx_t; \
    using MatrixNvNu_t = typename PhaseSpec::MatrixNvNu_t;   \
    using MatrixNuNv_t = typename PhaseSpec::MatrixNuNv_t;

#define GALILEO_PHASE_SPEC_ARRAY_TYPES_TYPEDEF(PhaseSpec)              \
    using VarScalarArray_t = typename PhaseSpec::VarScalarArray_t;     \
    using VectorNdxArray_t = typename PhaseSpec::VectorNdxArray_t;     \
    using VectorNxArray_t = typename PhaseSpec::VectorNxArray_t;       \
    using VectorNuArray_t = typename PhaseSpec::VectorNuArray_t;       \
    using MatrixNdxArray_t = typename PhaseSpec::MatrixNdxArray_t;     \
    using MatrixNdxNwArray_t = typename PhaseSpec::MatrixNdxNwArray_t; \
    using VectorNwArray_t = typename PhaseSpec::VectorNwArray_t;       \
    using MatrixNuArray_t = typename PhaseSpec::MatrixNuArray_t;       \
    using MatrixNdxNuArray_t = typename PhaseSpec::MatrixNdxNuArray_t; \
    using MatrixNuNwArray_t = typename PhaseSpec::MatrixNuNwArray_t;   \
    using MatrixNwArray_t = typename PhaseSpec::MatrixNwArray_t;

#define GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PhaseSpec)          \
    GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec)                \
    GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PhaseSpec)             \
    GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec)           \
    GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec)          \
    GALILEO_PHASE_SPEC_CONTROL_PARAM_TYPES_TYPEDEF(PhaseSpec) \
    GALILEO_PHASE_SPEC_SEGMENT_TYPES_TYPEDEF(PhaseSpec)       \
    GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec)         \
    GALILEO_PHASE_SPEC_ARRAY_TYPES_TYPEDEF(PhaseSpec)

namespace galileo
{

    namespace predictive
    {

        // // Declared elsewhere:
        // enum class NodeType
        // {
        //     FreeFwd = 0,    // Free forward dynamics node
        //     FreeInv = 1,    // Free inverse dynamics node
        //     ContactFwd = 2, // Contact forward dynamics node
        //     ContactInv = 3  // Contact inverse dynamics node
        // };

        // enum class SegmentType
        // {
        //     ERK = 0, // Explicit Runge-Kutta
        //     IRK = 1, // Implicit Runge-Kutta
        //     LIRK = 2 // Lifted Implicit Runge-Kutta
        // };

        /* ---------------------------------------------------------------- */
        /* Fully specifies the types and constants used in a phase. */
        /* ---------------------------------------------------------------- */
        template <typename BasicSpec,
                  template <typename> class ConstraintManagerTpl,
                  template <typename> class CostManagerTpl,
                  template <typename> class NodeTpl,
                  template <typename> class ControlParamTpl,
                  template <typename> class SegmentTpl,
                  template <typename> class PhaseTpl>
        struct PhaseSpecTpl
        {
            using BS = BasicSpec;
            using PS = PhaseSpecTpl<BS, ConstraintManagerTpl, CostManagerTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl>;

            /* ---------------------------------------------------------------- */
            /* Meta template types */
            /* ---------------------------------------------------------------- */
            using RobotMeta_t = typename BS::RobotMeta_t;
            using RobotModel_t = typename BS::RobotModel_t;
            using RobotData_t = typename BS::RobotData_t;

            using State_t = typename BS::State_t;

            using ActuationMeta_t = typename BS::ActuationMeta_t;
            using ActuationModel_t = typename BS::ActuationModel_t;
            using ActuationData_t = typename BS::ActuationData_t;

            using ConstraintManagerMeta_t = ConstraintManagerTpl<PS>;
            using ConstraintCollection_t = typename ConstraintManagerMeta_t::Collection;
            using ConstraintModelManager_t = typename ConstraintManagerMeta_t::Model;
            using ConstraintDataManager_t = typename ConstraintManagerMeta_t::Data;

            using CostManagerMeta_t = CostManagerTpl<PS>;
            using CostCollection_t = typename CostManagerMeta_t::Collection;
            using CostModelManager_t = typename CostManagerMeta_t::Model;
            using CostDataManager_t = typename CostManagerMeta_t::Data;

            using NodeMeta_t = NodeTpl<PS>;
            using NodeModel_t = typename NodeMeta_t::Model;
            using NodeData_t = typename NodeMeta_t::Data;
            using NodeDataVector_t = std::vector<NodeData_t>;

            using ControlParamMeta_t = ControlParamTpl<PS>;
            using ControlParamModel_t = typename ControlParamMeta_t::Model;
            using ControlParamData_t = typename ControlParamMeta_t::Data;
            using ControlParamDataVector_t = std::vector<ControlParamData_t>;

            using SegmentMeta_t = SegmentTpl<PS>;
            using SegmentModel_t = typename SegmentMeta_t::Model;
            using SegmentData_t = typename SegmentMeta_t::Data;
            using SegmentDataVector_t = std::vector<SegmentData_t>;

            using PhaseMeta_t = PhaseTpl<PS>;
            using PhaseModel_t = typename PhaseMeta_t::Model;
            using PhaseData_t = typename PhaseMeta_t::Data;
            using PhaseDataVector_t = std::vector<PhaseData_t>;

            /* ---------------------------------------------------------------- */
            /* Scalar types and Eigen Matrix storage order */
            /* ---------------------------------------------------------------- */
            using VarScalar = typename BS::VarScalar;   // Scalar type for variables (for AD)
            using NumScalar = typename BS::NumScalar;   // Scalar type for numerics (i.e., bounds, times, etc.)
            static constexpr int Options = BS::Options; // Eigen storage order

            /* ---------------------------------------------------------------- */
            /* Compile-time constants */
            /* ---------------------------------------------------------------- */
            static constexpr int NX = BS::NX;   // State dimension
            static constexpr int NDX = BS::NDX; // State tangent space dimension
            static constexpr int NQ = BS::NQ;   // Dimension of generalized coordinates
            static constexpr int NV = BS::NV;   // Dimension of generalized velocities

            /* ---------------------------------------------------------------- */
            /* Dependent compile-time constants */
            /* ---------------------------------------------------------------- */
            static constexpr int NU = traits<NodeMeta_t>::NU;                              // Control dimension
            static constexpr NodeType NodeType = traits<NodeMeta_t>::NodeType;             // Type of node: FreeFwd, FreeInv, ContactFwd, ContactInv
            static constexpr int NOrder = traits<ControlParamMeta_t>::NOrder;              // Order of the control parameterization per segment
            static constexpr int NW = traits<ControlParamMeta_t>::NW;                      // Number of control parameters
            static constexpr int NStages = traits<SegmentMeta_t>::NStages;                 // Number of Runge-Kutta stages per segment
            static constexpr SegmentType SegmentType = traits<SegmentMeta_t>::SegmentType; // Type of segment: ERK, IRK, LIRK

            /* ---------------------------------------------------------------- */
            /* Node type definitions */
            /* ---------------------------------------------------------------- */
            // Dynamics
            using XAcc_t = typename BS::VectorNv_t;     // System acceleration
            using XAccx_t = typename BS::MatrixNvNdx_t; // Jacobian of system acceleration w.r.t. state
            using XAccu_t = typename BS::MatrixNvNu_t;  // Jacobian of system acceleration w.r.t. control

            // Cost (CostManager holds an Eigen map to these, which are stored in NodeData)
            using L_t = VarScalar;                    // Cost scalar
            using Lx_t = typename BS::VectorNdx_t;    // Jacobian of cost w.r.t. state
            using Lu_t = typename BS::VectorNu_t;     // Jacobian of cost w.r.t. control
            using Lxx_t = typename BS::MatrixNdx_t;   // Hessian of cost w.r.t. state
            using Lxu_t = typename BS::MatrixNdxNu_t; // Hessian of cost w.r.t. state and control
            using Luu_t = typename BS::MatrixNu_t;    // Hessian of cost w.r.t. control

            // Equality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using H_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;          // Equality constraint vector
            using Hx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NDX, Options>;       // Jacobian of equality constraints w.r.t. state
            using Hu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>;        // Jacobian of equality constraints w.r.t. control
            using H_Equality_t = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1, Options>; // Equality constraint vector

            // Inequality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using G_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;       // Inequality constraint vector
            using Gx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NDX, Options>;    // Jacobian of inequality constraints w.r.t. state
            using Gu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>;     // Jacobian of inequality constraints w.r.t. control
            using G_Bound_t = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1, Options>; // Bounds on inequality constraints

            /* ---------------------------------------------------------------- */
            /* Control parameter type definitions */
            /* ---------------------------------------------------------------- */
            using U_t = typename BS::VectorNu_t;                    // Control vector
            using W_t = Eigen::Matrix<VarScalar, NW, 1, Options>;   // Control parameter vector
            using Uw_t = Eigen::Matrix<VarScalar, NU, NW, Options>; // Jacobian of control w.r.t. control parameters

            /* ---------------------------------------------------------------- */
            /* Segment type definitions */
            /* ---------------------------------------------------------------- */
            using Timings_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using Quadrature_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using StageCoefficients_t = Eigen::Matrix<NumScalar, NStages, NStages, Options>;

            // Dynamics
            using XNext_t = typename BS::VectorNx_t;                 // Evolution state
            using Fx_t = typename BS::MatrixNdx_t;                   // Jacobian of dynamics w.r.t. state
            using Fw_t = Eigen::Matrix<VarScalar, NDX, NW, Options>; // Jacobian of dynamics w.r.t. control parameters

            // Cost derivatives
            using Lw_t = Eigen::Matrix<VarScalar, NW, 1, Options>;    // Jacobian of cost w.r.t. control parameters
            using Lxw_t = Eigen::Matrix<VarScalar, NDX, NW, Options>; // Hessian of cost w.r.t. state and control parameters
            using Lww_t = Eigen::Matrix<VarScalar, NW, NW, Options>;  // Hessian of cost w.r.t. control parameters

            // Segment equality constraint derivatives
            using Hw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of equality constraints w.r.t. the control parameters

            // Segment inequality constraint derivatives
            using Gw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of inequality constraints w.r.t. the control parameters

            /* ---------------------------------------------------------------- */
            /* An assortment of Eigen types (primarily for use in Segments) */
            /* ---------------------------------------------------------------- */
            using VectorNx_t = typename BS::VectorNx_t;
            using VectorNu_t = typename BS::VectorNu_t;
            using VectorNdx_t = typename BS::VectorNdx_t;
            using VectorNq_t = typename BS::VectorNq_t;
            using VectorNv_t = typename BS::VectorNv_t;

            using MatrixNx_t = typename BS::MatrixNx_t;
            using MatrixNu_t = typename BS::MatrixNu_t;
            using MatrixNdx_t = typename BS::MatrixNdx_t;
            using MatrixNv_t = typename BS::MatrixNv_t;
            using MatrixNvNw_t = Eigen::Matrix<VarScalar, NV, NW, Options>;
            using MatrixNvNdx_t = typename BS::MatrixNvNdx_t;
            using MatrixNvNu_t = typename BS::MatrixNvNu_t;
            using MatrixNuNv_t = typename BS::MatrixNuNv_t;

            using VarScalarArray_t = std::array<VarScalar, NStages>;

            using VectorNdxArray_t = std::array<typename BS::VectorNdx_t, NStages>;
            using VectorNxArray_t = std::array<typename BS::VectorNx_t, NStages>;
            using VectorNuArray_t = std::array<typename BS::VectorNu_t, NStages>;
            using MatrixNdxArray_t = std::array<typename BS::MatrixNdx_t, NStages>;
            using MatrixNdxNwArray_t = std::array<Eigen::Matrix<VarScalar, NDX, NW, Options>, NStages>;
            using VectorNwArray_t = std::array<Eigen::Matrix<VarScalar, NW, 1, Options>, NStages>;
            using MatrixNuArray_t = std::array<typename BS::MatrixNu_t, NStages>;
            using MatrixNdxNuArray_t = std::array<typename BS::MatrixNdxNu_t, NStages>;
            using MatrixNuNwArray_t = std::array<Eigen::Matrix<VarScalar, NU, NW, Options>, NStages>;
            using MatrixNwArray_t = std::array<Eigen::Matrix<VarScalar, NW, NW, Options>, NStages>;
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_spec_hpp__
