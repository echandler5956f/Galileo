#ifndef __galileo_predictive_phases_phase_spec_hpp__
#define __galileo_predictive_phases_phase_spec_hpp__

#include "galileo/core/system-spec.hpp"
#include "galileo/predictive/phases/fwd.hpp"

#include <array>
#include <vector>

#define GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec) \
    using ConstraintManagerMeta_t = typename PhaseSpec::ConstraintManagerMeta_t; \
    using ConstraintCollection_t = typename PhaseSpec::ConstraintCollection_t; \
    using ConstraintModelManager_t = typename PhaseSpec::ConstraintModelManager_t; \
    using ConstraintDataManager_t = typename PhaseSpec::ConstraintDataManager_t; \
    using CostManagerMeta_t = typename PhaseSpec::CostManagerMeta_t; \
    using CostCollection_t = typename PhaseSpec::CostCollection_t; \
    using CostModelManager_t = typename PhaseSpec::CostModelManager_t; \
    using CostDataManager_t = typename PhaseSpec::CostDataManager_t; \
    using NodeMeta_t = typename PhaseSpec::NodeMeta_t; \
    using NodeModel_t = typename PhaseSpec::NodeModel_t; \
    using NodeData_t = typename PhaseSpec::NodeData_t; \
    using NodeModelVector_t = typename PhaseSpec::NodeModelVector_t; \
    using NodeDataVector_t = typename PhaseSpec::NodeDataVector_t; \
    using ControlParamMeta_t = typename PhaseSpec::ControlParamMeta_t; \
    using ControlParamModel_t = typename PhaseSpec::ControlParamModel_t; \
    using ControlParamData_t = typename PhaseSpec::ControlParamData_t; \
    using ControlParamModelVector_t = typename PhaseSpec::ControlParamModelVector_t; \
    using ControlParamDataVector_t = typename PhaseSpec::ControlParamDataVector_t; \
    using SegmentMeta_t = typename PhaseSpec::SegmentMeta_t; \
    using SegmentModel_t = typename PhaseSpec::SegmentModel_t; \
    using SegmentData_t = typename PhaseSpec::SegmentData_t; \
    using SegmentModelVector_t = typename PhaseSpec::SegmentModelVector_t; \
    using SegmentDataVector_t = typename PhaseSpec::SegmentDataVector_t; \
    using PhaseMeta_t = typename PhaseSpec::PhaseMeta_t; \
    using PhaseModel_t = typename PhaseSpec::PhaseModel_t; \
    using PhaseData_t = typename PhaseSpec::PhaseData_t;

#define GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PhaseSpec) GALILEO_SYSTEM_SPEC_SCALARS_TYPEDEF(PhaseSpec::SS);

#define GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec) \
    GALILEO_SYSTEM_SPEC_CONSTANTS_TYPEDEF(PhaseSpec::SS); \
    static constexpr int NU = PhaseSpec::NU; \
    static constexpr int NOrder = PhaseSpec::NOrder; \
    static constexpr int NW = PhaseSpec::NW; \
    static constexpr int NStages = PhaseSpec::NStages; \
    using DimNU_t = typename PhaseSpec::DimNU_t; \
    using DimNOrder_t = typename PhaseSpec::DimNOrder_t; \
    using DimNW_t = typename PhaseSpec::DimNW_t; \
    using DimNStages_t = typename PhaseSpec::DimNStages_t;

#define GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec) \
    using XAcc_t = typename PhaseSpec::XAcc_t; \
    using XAccx_t = typename PhaseSpec::XAccx_t; \
    using XAccu_t = typename PhaseSpec::XAccu_t; \
    using L_t = typename PhaseSpec::L_t; \
    using Lx_t = typename PhaseSpec::Lx_t; \
    using Lu_t = typename PhaseSpec::Lu_t; \
    using Lxx_t = typename PhaseSpec::Lxx_t; \
    using Lxu_t = typename PhaseSpec::Lxu_t; \
    using Luu_t = typename PhaseSpec::Luu_t; \
    using H_t = typename PhaseSpec::H_t; \
    using Hx_t = typename PhaseSpec::Hx_t; \
    using Hu_t = typename PhaseSpec::Hu_t; \
    using G_t = typename PhaseSpec::G_t; \
    using Gx_t = typename PhaseSpec::Gx_t; \
    using Gu_t = typename PhaseSpec::Gu_t; \
    using G_Bound_t = typename PhaseSpec::G_Bound_t;

#define GALILEO_PHASE_SPEC_CONTROL_PARAM_TYPES_TYPEDEF(PhaseSpec) \
    using U_t = typename PhaseSpec::U_t; \
    using W_t = typename PhaseSpec::W_t; \
    using Uw_t = typename PhaseSpec::Uw_t;

#define GALILEO_PHASE_SPEC_SEGMENT_TYPES_TYPEDEF(PhaseSpec) \
    using Timings_t = typename PhaseSpec::Timings_t; \
    using Quadrature_t = typename PhaseSpec::Quadrature_t; \
    using StageCoefficients_t = typename PhaseSpec::StageCoefficients_t; \
    using XNext_t = typename PhaseSpec::XNext_t; \
    using XNextx_t = typename PhaseSpec::XNextx_t; \
    using XNextw_t = typename PhaseSpec::XNextw_t; \
    using Lw_t = typename PhaseSpec::Lw_t; \
    using Lxw_t = typename PhaseSpec::Lxw_t; \
    using Lww_t = typename PhaseSpec::Lww_t; \
    using Hw_t = typename PhaseSpec::Hw_t; \
    using Gw_t = typename PhaseSpec::Gw_t;

#define GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec) \
    GALILEO_SYSTEM_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec::SS); \
    using VectorNu_t = typename PhaseSpec::VectorNu_t; \
    using VectorNw_t = typename PhaseSpec::VectorNw_t; \
    using MatrixNu_t = typename PhaseSpec::MatrixNu_t; \
    using MatrixNw_t = typename PhaseSpec::MatrixNw_t; \
    using MatrixNvNw_t = typename PhaseSpec::MatrixNvNw_t; \
    using MatrixNvNu_t = typename PhaseSpec::MatrixNvNu_t; \
    using MatrixNuNv_t = typename PhaseSpec::MatrixNuNv_t; \
    using MatrixNuNw_t = typename PhaseSpec::MatrixNuNw_t; \
    using MatrixNuaNu_t = typename PhaseSpec::MatrixNuaNu_t; \
    using MatrixNdxNw_t = typename PhaseSpec::MatrixNdxNw_t;

#define GALILEO_PHASE_SPEC_ARRAY_TYPES_TYPEDEF(PhaseSpec) \
    using VarScalarArray_t = typename PhaseSpec::VarScalarArray_t; \
    using VectorNdxArray_t = typename PhaseSpec::VectorNdxArray_t; \
    using VectorNxArray_t = typename PhaseSpec::VectorNxArray_t; \
    using VectorNuArray_t = typename PhaseSpec::VectorNuArray_t; \
    using MatrixNdxArray_t = typename PhaseSpec::MatrixNdxArray_t; \
    using MatrixNdxNwArray_t = typename PhaseSpec::MatrixNdxNwArray_t; \
    using VectorNwArray_t = typename PhaseSpec::VectorNwArray_t; \
    using MatrixNuArray_t = typename PhaseSpec::MatrixNuArray_t; \
    using MatrixNdxNuArray_t = typename PhaseSpec::MatrixNdxNuArray_t; \
    using MatrixNuNwArray_t = typename PhaseSpec::MatrixNuNwArray_t; \
    using MatrixNwArray_t = typename PhaseSpec::MatrixNwArray_t;

#define GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PhaseSpec) \
    GALILEO_PHASE_SPEC_META_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_SCALARS_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_CONSTANTS_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_NODE_TYPES_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_CONTROL_PARAM_TYPES_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_SEGMENT_TYPES_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_EIGEN_TYPES_TYPEDEF(PhaseSpec); \
    GALILEO_PHASE_SPEC_ARRAY_TYPES_TYPEDEF(PhaseSpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Fully specifies the types and constants used in a phase. */
    /* ---------------------------------------------------------------- */
    template <typename DomainSystemSpec,
              template <typename> class ConstraintManagerTpl,
              template <typename> class CostManagerTpl,
              template <typename> class NodeTpl,
              template <typename> class ControlParamTpl,
              template <typename> class SegmentTpl,
              template <typename> class PhaseTpl>
    struct PhaseSpecTpl : public DomainSystemSpec
    {
        using PS = PhaseSpecTpl<DomainSystemSpec,
                                ConstraintManagerTpl,
                                CostManagerTpl,
                                NodeTpl,
                                ControlParamTpl,
                                SegmentTpl,
                                PhaseTpl>;
        using SS = DomainSystemSpec;
        using Base = DomainSystemSpec;
        using BS = typename SS::BS;

        /* ---------------------------------------------------------------- */
        /* Import the system spec types and constants */
        /* ---------------------------------------------------------------- */
        GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(Base);

        using State_t = typename Base::State_t;

        using ActuationMeta_t = typename Base::ActuationMeta_t;
        using ActuationModel_t = typename Base::ActuationModel_t;
        using ActuationData_t = typename Base::ActuationData_t;

        /* ---------------------------------------------------------------- */
        /* Meta template types */
        /* ---------------------------------------------------------------- */
        using NodeMeta_t = NodeTpl<PS>;
        using NodeModel_t = typename traits<NodeMeta_t>::Model_t;
        using NodeData_t = typename traits<NodeMeta_t>::Data_t;

        using NodeModelVector_t = std::vector<NodeModel_t>;
        using NodeDataVector_t = std::vector<NodeData_t>;

        using ControlParamMeta_t = ControlParamTpl<PS>;
        using ControlParamModel_t = typename traits<ControlParamMeta_t>::Model_t;
        using ControlParamData_t = typename traits<ControlParamMeta_t>::Data_t;

        using ControlParamModelVector_t = std::vector<ControlParamModel_t>;
        using ControlParamDataVector_t = std::vector<ControlParamData_t>;

        using SegmentMeta_t = SegmentTpl<PS>;
        using SegmentModel_t = typename traits<SegmentMeta_t>::Model_t;
        using SegmentData_t = typename traits<SegmentMeta_t>::Data_t;

        using SegmentModelVector_t = std::vector<SegmentModel_t>;
        using SegmentDataVector_t = std::vector<SegmentData_t>;

        /* ---------------------------------------------------------------- */
        /* Dependent dimension types */
        /* ---------------------------------------------------------------- */

        // ANYTHING THAT USES THE DEPENDENT CONSTANTS MUST BE DEFINED
        // *AFTER* THESE DEPENDENT CONSTANTS ARE DEFINED

        using DimNU_t = typename traits<NodeMeta_t>::DimNU_t;
        using DimNOrder_t = typename traits<ControlParamMeta_t>::DimNOrder_t;
        using DimNW_t = decltype(DimNU_t{} * DimNOrder_t{});
        using DimNStages_t = typename traits<SegmentMeta_t>::DimNStages_t;

        /* NOTE: NU is calculated differently depending on the node type */
        // FreeFwd: NU = NUa
        // FreeInv: NU = NV
        // ContactFwd: NU = NUa
        // ContactInv: NU = NV + NContacts
        // We do this calculation in the traits specialization for each derived node type.

        /* ---------------------------------------------------------------- */
        /* Compile-time constants */
        /* ---------------------------------------------------------------- */
        static constexpr int NU = DimNU_t::Value;           // Control dimension
        static constexpr int NOrder = DimNOrder_t::Value;   // Order of the control parameterization per segment
        static constexpr int NW = DimNW_t::Value;           // Number of control parameters
        static constexpr int NStages = DimNStages_t::Value; // Number of Runge-Kutta stages per segment

        using ConstraintManagerMeta_t = ConstraintManagerTpl<PS>;
        using ConstraintCollection_t = typename traits<ConstraintManagerMeta_t>::Collection_t;
        using ConstraintModelManager_t = typename traits<ConstraintManagerMeta_t>::ModelManager_t;
        using ConstraintDataManager_t = typename traits<ConstraintManagerMeta_t>::DataManager_t;

        using CostManagerMeta_t = CostManagerTpl<PS>;
        using CostCollection_t = typename traits<CostManagerMeta_t>::Collection_t;
        using CostModelManager_t = typename traits<CostManagerMeta_t>::ModelManager_t;
        using CostDataManager_t = typename traits<CostManagerMeta_t>::DataManager_t;

        /* ---------------------------------------------------------------- */
        /* An assortment of Eigen types (primarily for use in Segments) */
        /* ---------------------------------------------------------------- */
        using VectorNu_t = Eigen::GMatrix<VarScalar, NU, 1, Options>;
        using VectorNw_t = Eigen::GMatrix<VarScalar, NW, 1, Options>;
        using MatrixNu_t = Eigen::GMatrix<VarScalar, NU, NU, Options>;
        using MatrixNw_t = Eigen::GMatrix<VarScalar, NW, NW, Options>;

        using MatrixNvNw_t = Eigen::GMatrix<VarScalar, NV, NW, Options>;
        using MatrixNvNu_t = Eigen::GMatrix<VarScalar, NV, NU, Options>;

        using MatrixNuNv_t = Eigen::GMatrix<VarScalar, NU, NV, Options>;
        using MatrixNuNw_t = Eigen::GMatrix<VarScalar, NU, NW, Options>;

        using MatrixNuaNu_t = Eigen::GMatrix<VarScalar, NUa, NU, Options>;

        using MatrixNdxNu_t = Eigen::GMatrix<VarScalar, NDX, NU, Options>;
        using MatrixNdxNw_t = Eigen::GMatrix<VarScalar, NDX, NW, Options>;

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
        using H_t = Eigen::GMatrix<VarScalar, Eigen::Dynamic, 1, Options>; // Equality constraint vector
        using Hx_t =
            Eigen::GMatrix<VarScalar, Eigen::Dynamic, NDX, Options>; // Jacobian of equality constraints w.r.t. state
        using Hu_t =
            Eigen::GMatrix<VarScalar, Eigen::Dynamic, NU, Options>; // Jacobian of equality constraints w.r.t. control

        // Inequality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
        using G_t = Eigen::GMatrix<VarScalar, Eigen::Dynamic, 1, Options>; // Inequality constraint vector
        using Gx_t =
            Eigen::GMatrix<VarScalar, Eigen::Dynamic, NDX, Options>; // Jacobian of inequality constraints w.r.t. state
        using Gu_t =
            Eigen::GMatrix<VarScalar, Eigen::Dynamic, NU, Options>; // Jacobian of inequality constraints w.r.t. control
        using G_Bound_t = Eigen::GMatrix<NumScalar, Eigen::Dynamic, 1, Options>; // Bounds on inequality constraints

        /* ---------------------------------------------------------------- */
        /* Control parameter type definitions */
        /* ---------------------------------------------------------------- */
        using U_t = VectorNu_t;    // Control vector
        using W_t = VectorNw_t;    // Control parameter vector
        using Uw_t = MatrixNuNw_t; // Jacobian of control w.r.t. control parameters

        /* ---------------------------------------------------------------- */
        /* Segment type definitions */
        /* ---------------------------------------------------------------- */
        using Timings_t = Eigen::GMatrix<NumScalar, NStages, 1, Options>;
        using Quadrature_t = Eigen::GMatrix<NumScalar, NStages, 1, Options>;
        using StageCoefficients_t = Eigen::GMatrix<NumScalar, NStages, NStages, Options>;

        // Dynamics
        using XNext_t = VectorNx_t;     // Evolution state
        using XNextx_t = MatrixNdx_t;   // Jacobian of dynamics w.r.t. state
        using XNextw_t = MatrixNdxNw_t; // Jacobian of dynamics w.r.t. control parameters

        // Cost derivatives
        using Lw_t = VectorNw_t;     // Jacobian of cost w.r.t. control parameters
        using Lxw_t = MatrixNdxNw_t; // Hessian of cost w.r.t. state and control parameters
        using Lww_t = MatrixNw_t;    // Hessian of cost w.r.t. control parameters

        // Segment equality constraint derivatives
        using Hw_t = Eigen::GMatrix<VarScalar,
                                    Eigen::Dynamic,
                                    NW,
                                    Options>; // Jacobian of equality constraints w.r.t. the control parameters

        // Segment inequality constraint derivatives
        using Gw_t = Eigen::GMatrix<VarScalar,
                                    Eigen::Dynamic,
                                    NW,
                                    Options>; // Jacobian of inequality constraints w.r.t. the control parameters

        using PhaseMeta_t = PhaseTpl<PS>;
        using PhaseModel_t = typename traits<PhaseMeta_t>::Model_t;
        using PhaseData_t = typename traits<PhaseMeta_t>::Data_t;

        /* ---------------------------------------------------------------- */
        /* Additional dimension storage */
        /* ---------------------------------------------------------------- */
        using Base::nq_dim_;
        using Base::nv_dim_;
        using Base::nx_dim_;
        using Base::ndx_dim_;
        using Base::nua_dim_;
        using Base::nqb_dim_;
        using Base::nqj_dim_;
        using Base::nvb_dim_;
        using Base::nvj_dim_;
        using Base::nrotors_dim_;

        DimNU_t nu_dim_;
        DimNOrder_t norder_dim_;
        DimNW_t nw_dim_;
        DimNStages_t nstages_dim_;

        /* ---------------------------------------------------------------- */
        /* Constructors */
        /* ---------------------------------------------------------------- */
        PhaseSpecTpl(const DomainSystemSpec &ss, int NU_ = NU, int NOrder_ = NOrder, int NStages_ = NStages)
            : Base(ss), nu_dim_{NU_}, norder_dim_{NOrder_}, nw_dim_(nu_dim_ * norder_dim_), nstages_dim_{NStages_}
        {
        }

        /* ---------------------------------------------------------------- */
        /* Accessors for the PhaseSpec dimensions */
        /* ---------------------------------------------------------------- */
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

        int get_nu() const { return nu_dim_.value(); }
        const DimNU_t &get_nu_dim() const { return nu_dim_; }
        int get_norder() const { return norder_dim_.value(); }
        const DimNOrder_t &get_norder_dim() const { return norder_dim_; }
        int get_nw() const { return nw_dim_.value(); }
        const DimNW_t &get_nw_dim() const { return nw_dim_; }
        int get_nstages() const { return nstages_dim_.value(); }
        const DimNStages_t &get_nstages_dim() const { return nstages_dim_; }

        inline bool is_valid_spec() const
        {
            bool valid_base = Base::is_valid_spec();
            bool valid_nu = (get_nu() >= 0);
            bool valid_norder = (get_norder() >= 0);
            bool valid_nw = (get_nw() == get_nu() * get_norder());
            bool valid_nstages = (get_nstages() > 0); // Zero stages is not allowed
            return valid_base && valid_nu && valid_norder && valid_nw && valid_nstages;
        }

        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            Base::display(os, indent);

            os << indent << "Phase-Specific Dimensions: {\n";
            os << indent << "  Control: {\n";
            os << indent << "    NU (Control dimension):         " << nu_dim_ << "\n";
            os << indent << "    NOrder (Control order):         " << norder_dim_ << "\n";
            os << indent << "    NW (Control parameters):        " << nw_dim_ << "\n";
            os << indent << "  }\n";
            os << indent << "  Integration: {\n";
            os << indent << "    NStages (Runge-Kutta stages):   " << nstages_dim_ << "\n";
            os << indent << "  }\n";
            os << indent << "}\n";
            os << indent << "Spec configuration is " << (is_valid_spec() ? "VALID" : "INVALID") << "\n";
        }

        // Stream output operator
        friend std::ostream &operator<<(std::ostream &os, const PhaseSpecTpl &ps)
        {
            os << "PhaseSpec: {\n";
            ps.display(os, "  ");
            os << "}";
            return os;
        }
    };

} // namespace galileo

#endif // __galileo_predictive_phases_phase_spec_hpp__
