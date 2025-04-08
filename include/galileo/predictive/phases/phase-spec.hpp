#ifndef __galileo_predictive_phases_phase_spec_hpp__
#define __galileo_predictive_phases_phase_spec_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include "galileo/core/basic-spec.hpp"

#include <vector>
#include <array>

#define GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec)                                 \
    GALILEO_BASIC_SPEC_META_TYPEDEF(PhaseSpec::BS);                                \
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
    GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(PhaseSpec::BS);

#define GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec)  \
    GALILEO_BASIC_SPEC_CONSTANTS_TYPEDEF(PhaseSpec::BS); \
    static constexpr int NU = PhaseSpec::NU;             \
    static constexpr int NOrder = PhaseSpec::NOrder;     \
    static constexpr int NW = PhaseSpec::NW;             \
    static constexpr int NStages = PhaseSpec::NStages;

#define GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec) \
    using XAcc_t = typename PhaseSpec::XAcc_t;           \
    using XAccx_t = typename PhaseSpec::XAccx_t;         \
    using XAccu_t = typename PhaseSpec::XAccu_t;         \
    using L_t = typename PhaseSpec::L_t;                 \
    using Lx_t = typename PhaseSpec::Lx_t;               \
    using Lu_t = typename PhaseSpec::Lu_t;               \
    using Lxx_t = typename PhaseSpec::Lxx_t;             \
    using Lxu_t = typename PhaseSpec::Lxu_t;             \
    using Luu_t = typename PhaseSpec::Luu_t;             \
    using H_t = typename PhaseSpec::H_t;                 \
    using Hx_t = typename PhaseSpec::Hx_t;               \
    using Hu_t = typename PhaseSpec::Hu_t;               \
    using G_t = typename PhaseSpec::G_t;                 \
    using Gx_t = typename PhaseSpec::Gx_t;               \
    using Gu_t = typename PhaseSpec::Gu_t;               \
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
    using XNextx_t = typename PhaseSpec::XNextx_t;                       \
    using XNextw_t = typename PhaseSpec::XNextw_t;                       \
    using Lw_t = typename PhaseSpec::Lw_t;                               \
    using Lxw_t = typename PhaseSpec::Lxw_t;                             \
    using Lww_t = typename PhaseSpec::Lww_t;                             \
    using Hw_t = typename PhaseSpec::Hw_t;                               \
    using Gw_t = typename PhaseSpec::Gw_t;

#define GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec)    \
    GALILEO_BASIC_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec::BS);   \
    using VectorNu_t = typename PhaseSpec::VectorNu_t;       \
    using VectorNw_t = typename PhaseSpec::VectorNw_t;       \
    using MatrixNu_t = typename PhaseSpec::MatrixNu_t;       \
    using MatrixNw_t = typename PhaseSpec::MatrixNw_t;       \
    using MatrixNvNw_t = typename PhaseSpec::MatrixNvNw_t;   \
    using MatrixNvNu_t = typename PhaseSpec::MatrixNvNu_t;   \
    using MatrixNvNua_t = typename PhaseSpec::MatrixNvNua_t; \
    using MatrixNvNdx_t = typename PhaseSpec::MatrixNvNdx_t; \
    using MatrixNuNv_t = typename PhaseSpec::MatrixNuNv_t;   \
    using MatrixNuNw_t = typename PhaseSpec::MatrixNuNw_t;   \
    using MatrixNuaNu_t = typename PhaseSpec::MatrixNuaNu_t; \
    using MatrixNuaNv_t = typename PhaseSpec::MatrixNuaNv_t; \
    using MatrixNdxNw_t = typename PhaseSpec::MatrixNdxNw_t; \
    using MatrixNdxNua_t = typename PhaseSpec::MatrixNdxNua_t;

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

#define GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PhaseSpec)           \
    GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec);                \
    GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PhaseSpec);             \
    GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec);           \
    GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec);          \
    GALILEO_PHASE_SPEC_CONTROL_PARAM_TYPES_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_SEGMENT_TYPES_TYPEDEF(PhaseSpec);       \
    GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec);         \
    GALILEO_PHASE_SPEC_ARRAY_TYPES_TYPEDEF(PhaseSpec);

namespace galileo
{

    namespace predictive
    {

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

            // Import the basic spec types and constants
            GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

            /* ---------------------------------------------------------------- */
            /* Meta template types */
            /* ---------------------------------------------------------------- */
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
            /* Dependent compile-time constants */
            /* ---------------------------------------------------------------- */
            static constexpr int NU = traits<NodeMeta_t>::NU;                 // Control dimension
            static constexpr int NOrder = traits<ControlParamMeta_t>::NOrder; // Order of the control parameterization per segment
            static constexpr int NW = traits<ControlParamMeta_t>::NW;         // Number of control parameters
            static constexpr int NStages = traits<SegmentMeta_t>::NStages;    // Number of Runge-Kutta stages per segment

            /*NOTE: NU is calculated differently depending on the node type*/
            // FreeFwd: NU = NUa
            // FreeInv: NU = NV
            // ContactFwd: NU = NUa
            // ContactInv: NU = NV + NContacts
            // We do this calculation in the traits specialization for each derived node type.

            /* ---------------------------------------------------------------- */
            /* An assortment of Eigen types (primarily for use in Segments) */
            /* ---------------------------------------------------------------- */
            using VectorNu_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            using VectorNw_t = Eigen::Matrix<VarScalar, NW, 1, Options>;
            using MatrixNu_t = Eigen::Matrix<VarScalar, NU, NU, Options>;
            using MatrixNw_t = Eigen::Matrix<VarScalar, NW, NW, Options>;

            using MatrixNvNw_t = Eigen::Matrix<VarScalar, NV, NW, Options>;
            using MatrixNvNu_t = Eigen::Matrix<VarScalar, NV, NU, Options>;
            using MatrixNvNua_t = Eigen::Matrix<VarScalar, NV, NUa, Options>;
            using MatrixNvNdx_t = Eigen::Matrix<VarScalar, NV, NDX, Options>;

            using MatrixNuNv_t = Eigen::Matrix<VarScalar, NU, NV, Options>;
            using MatrixNuNw_t = Eigen::Matrix<VarScalar, NU, NW, Options>;

            using MatrixNuaNu_t = Eigen::Matrix<VarScalar, NUa, NU, Options>;
            using MatrixNuaNv_t = Eigen::Matrix<VarScalar, NUa, NV, Options>;

            using MatrixNdxNw_t = Eigen::Matrix<VarScalar, NDX, NW, Options>;
            using MatrixNdxNua_t = Eigen::Matrix<VarScalar, NDX, NUa, Options>;

            using VarScalarArray_t = std::array<VarScalar, NStages>;

            using VectorNdxArray_t = std::array<VectorNdx_t, NStages>;
            using VectorNxArray_t = std::array<VectorNx_t, NStages>;
            using VectorNuArray_t = std::array<VectorNu_t, NStages>;
            using MatrixNdxArray_t = std::array<MatrixNdx_t, NStages>;
            using MatrixNdxNwArray_t = std::array<MatrixNdxNw_t, NStages>;
            using VectorNwArray_t = std::array<VectorNw_t, NStages>;
            using MatrixNuArray_t = std::array<MatrixNu_t, NStages>;
            using MatrixNdxNuArray_t = std::array<MatrixNdxNu_t, NStages>;
            using MatrixNuNwArray_t = std::array<MatrixNuNw_t, NStages>;
            using MatrixNwArray_t = std::array<MatrixNw_t, NStages>;

            /* ---------------------------------------------------------------- */
            /* Node type definitions */
            /* ---------------------------------------------------------------- */
            // Dynamics
            using XAcc_t = VectorNv_t;     // System acceleration
            using XAccx_t = MatrixNvNdx_t; // Jacobian of system acceleration w.r.t. state
            using XAccu_t = MatrixNvNu_t;  // Jacobian of system acceleration w.r.t. control

            // Cost (CostManager holds an Eigen map to these, which are stored in NodeData)
            using L_t = VarScalar;       // Cost scalar
            using Lx_t = VectorNdx_t;    // Jacobian of cost w.r.t. state
            using Lu_t = VectorNu_t;     // Jacobian of cost w.r.t. control
            using Lxx_t = MatrixNdx_t;   // Hessian of cost w.r.t. state
            using Lxu_t = MatrixNdxNu_t; // Hessian of cost w.r.t. state and control
            using Luu_t = MatrixNu_t;    // Hessian of cost w.r.t. control

            // Equality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using H_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;    // Equality constraint vector
            using Hx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NDX, Options>; // Jacobian of equality constraints w.r.t. state
            using Hu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>;  // Jacobian of equality constraints w.r.t. control

            // Inequality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using G_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;       // Inequality constraint vector
            using Gx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NDX, Options>;    // Jacobian of inequality constraints w.r.t. state
            using Gu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>;     // Jacobian of inequality constraints w.r.t. control
            using G_Bound_t = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1, Options>; // Bounds on inequality constraints

            /* ---------------------------------------------------------------- */
            /* Control parameter type definitions */
            /* ---------------------------------------------------------------- */
            using U_t = VectorNu_t;    // Control vector
            using W_t = VectorNw_t;    // Control parameter vector
            using Uw_t = MatrixNuNw_t; // Jacobian of control w.r.t. control parameters

            /* ---------------------------------------------------------------- */
            /* Segment type definitions */
            /* ---------------------------------------------------------------- */
            using Timings_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using Quadrature_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using StageCoefficients_t = Eigen::Matrix<NumScalar, NStages, NStages, Options>;

            // Dynamics
            using XNext_t = VectorNx_t;     // Evolution state
            using XNextx_t = MatrixNdx_t;   // Jacobian of dynamics w.r.t. state
            using XNextw_t = MatrixNdxNw_t; // Jacobian of dynamics w.r.t. control parameters

            // Cost derivatives
            using Lw_t = VectorNw_t;     // Jacobian of cost w.r.t. control parameters
            using Lxw_t = MatrixNdxNw_t; // Hessian of cost w.r.t. state and control parameters
            using Lww_t = MatrixNw_t;    // Hessian of cost w.r.t. control parameters

            // Segment equality constraint derivatives
            using Hw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of equality constraints w.r.t. the control parameters

            // Segment inequality constraint derivatives
            using Gw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of inequality constraints w.r.t. the control parameters
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_spec_hpp__
