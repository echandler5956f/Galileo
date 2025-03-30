#ifndef __galileo_predictive_phases_phase_spec_hpp__
#define __galileo_predictive_phases_phase_spec_hpp__

#include "galileo/predictive/phases/fwd.hpp"

#include <vector>
#include <array>

namespace galileo
{

    namespace predictive
    {

        enum class SegmentType
        {
            ERK = 0, // Explicit Runge-Kutta
            IRK = 1, // Implicit Runge-Kutta
            LIRK = 2 // Lifted Implicit Runge-Kutta
        };

        /* ---------------------------------------------------------------- */
        /* Fully specifies the types and constants used in a phase. */
        /* ---------------------------------------------------------------- */
        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX,
                  int _NQ,
                  int _NV,
                  int _NStages,
                  int _NOrder,
                  template <typename> class StateTpl,
                  template <typename> class ActuationTpl,
                  template <typename> class ConstraintManagerTpl,
                  template <typename> class CostManagerTpl,
                  template <typename> class NodeTpl,
                  template <typename> class ControlParamTpl,
                  template <typename> class SegmentTpl,
                  template <typename> class PhaseTpl,
                  SegmentType _SegmentType>
        struct PhaseSpecTpl
        {
            using PhaseSpec = PhaseSpecTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX, _NQ, _NV, _NStages, _NOrder, StateTpl, ActuationTpl, ConstraintManagerTpl, CostManagerTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl, _SegmentType>;

            /* ---------------------------------------------------------------- */
            /* Scalar types and Eigen Matrix storage order */
            /* ---------------------------------------------------------------- */
            using VarScalar = _VarScalar;            // Scalar type for variables (for AD)
            using NumScalar = _NumScalar;            // Scalar type for numerics (i.e., bounds, times, etc.)
            static constexpr int Options = _Options; // Eigen storage order

            /* ---------------------------------------------------------------- */
            /* Compile-time constants */
            /* ---------------------------------------------------------------- */
            static constexpr int NX = _NX;           // State dimension
            static constexpr int NU = _NU;           // Control dimension
            static constexpr int NDX = _NDX;         // State tangent space dimension
            static constexpr int NQ = _NQ;           // Dimension of generalized coordinates
            static constexpr int NV = _NV;           // Dimension of generalized velocities
            static constexpr int NStages = _NStages; // Number of Runge-Kutta stages per segment
            static constexpr int NOrder = _NOrder;   // Order of the control parameterization per segment
            static constexpr int NW = NU * NOrder;   // Number of control parameters

            static constexpr SegmentType SegmentType = _SegmentType;

            /* ---------------------------------------------------------------- */
            /* Node type definitions */
            /* ---------------------------------------------------------------- */
            // Dynamics
            using XAcc_t = Eigen::Matrix<VarScalar, NV, 1, Options>;    // System acceleration
            using XAccx_t = Eigen::Matrix<VarScalar, NV, NDX, Options>; // Jacobian of system acceleration w.r.t. state
            using XAccu_t = Eigen::Matrix<VarScalar, NV, NU, Options>;  // Jacobian of system acceleration w.r.t. control

            // Cost (CostManager holds an Eigen map to these, which are stored in NodeData)
            using L_t = VarScalar;                                     // Cost scalar
            using Lx_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;    // Jacobian of cost w.r.t. state
            using Lu_t = Eigen::Matrix<VarScalar, NU, 1, Options>;     // Jacobian of cost w.r.t. control
            using Lxx_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>; // Hessian of cost w.r.t. state
            using Lxu_t = Eigen::Matrix<VarScalar, NDX, NU, Options>;  // Hessian of cost w.r.t. state and control
            using Luu_t = Eigen::Matrix<VarScalar, NU, NU, Options>;   // Hessian of cost w.r.t. control

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
            using U_t = Eigen::Matrix<NumScalar, NU, 1, Options>;   // Control vector
            using W_t = Eigen::Matrix<NumScalar, NW, 1, Options>;   // Control parameter vector
            using Uw_t = Eigen::Matrix<NumScalar, NU, NW, Options>; // Jacobian of control w.r.t. control parameters

            /* ---------------------------------------------------------------- */
            /* Segment type definitions */
            /* ---------------------------------------------------------------- */
            using Timings_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using Quadrature_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using StageCoefficients_t = Eigen::Matrix<NumScalar, NStages, NStages, Options>;

            // Dynamics
            using XNext_t = Eigen::Matrix<VarScalar, NX, 1, Options>; // Evolution state
            using Fx_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>; // Jacobian of dynamics w.r.t. state
            using Fw_t = Eigen::Matrix<VarScalar, NDX, NW, Options>;  // Jacobian of dynamics w.r.t. control parameters

            // Cost derivatives
            using Lw_t = Eigen::Matrix<VarScalar, NU, 1, Options>;    // Jacobian of cost w.r.t. control parameters
            using Lxw_t = Eigen::Matrix<VarScalar, NDX, NW, Options>; // Hessian of cost w.r.t. state and control parameters
            using Lww_t = Eigen::Matrix<VarScalar, NW, NW, Options>;  // Hessian of cost w.r.t. control parameters

            // Segment equality constraint derivatives
            using Hw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of equality constraints w.r.t. the control parameters

            // Segment inequality constraint derivatives
            using Gw_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NW, Options>; // Jacobian of inequality constraints w.r.t. the control parameters

            /* ---------------------------------------------------------------- */
            /* An assortment of Eigen types (primarily for use in Segments) */
            /* ---------------------------------------------------------------- */
            using VectorNx_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
            using VectorNu_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            /*1*/ using VectorNdx_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using VectorNq_t = Eigen::Matrix<VarScalar, NQ, 1, Options>;
            using VectorNv_t = Eigen::Matrix<VarScalar, NV, 1, Options>;

            using MatrixNx_t = Eigen::Matrix<VarScalar, NX, NX, Options>;
            using MatrixNu_t = Eigen::Matrix<VarScalar, NU, NU, Options>;
            using MatrixNdx_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
            using MatrixNv_t = Eigen::Matrix<VarScalar, NV, NV, Options>;
            using MatrixNvNw_t = Eigen::Matrix<VarScalar, NV, NW, Options>;
            using MatrixNvNdx_t = Eigen::Matrix<VarScalar, NV, NDX, Options>;
            using MatrixNvNu_t = Eigen::Matrix<VarScalar, NV, NU, Options>;
            using MatrixNuNv_t = Eigen::Matrix<VarScalar, NU, NV, Options>;

            using VarScalarArray_t = std::array<VarScalar, NStages>;

            /*2*/ using VectorNdxArray_t = std::array<Eigen::Matrix<VarScalar, NDX, 1, Options>, NStages>;
            /*3*/ using VectorNxArray_t = std::array<Eigen::Matrix<VarScalar, NX, 1, Options>, NStages>;
            /*4*/ using VectorNuArray_t = std::array<Eigen::Matrix<VarScalar, NU, 1, Options>, NStages>;
            /*5*/ using MatrixNdxArray_t = std::array<Eigen::Matrix<VarScalar, NDX, NDX, Options>, NStages>;
            /*6*/ using MatrixNdxNwArray_t = std::array<Eigen::Matrix<VarScalar, NDX, NW, Options>, NStages>;
            /*7*/ using VectorNwArray_t = std::array<Eigen::Matrix<VarScalar, NW, 1, Options>, NStages>;
            /*8*/ using MatrixNuArray_t = std::array<Eigen::Matrix<VarScalar, NU, NU, Options>, NStages>;
            /*9*/ using MatrixNdxNuArray_t = std::array<Eigen::Matrix<VarScalar, NDX, NU, Options>, NStages>;
            /*10*/ using MatrixNuNwArray_t = std::array<Eigen::Matrix<VarScalar, NU, NW, Options>, NStages>;
            /*11*/ using MatrixNwArray_t = std::array<Eigen::Matrix<VarScalar, NW, NW, Options>, NStages>;

            // 1  (ndx, 1)           // dx
            // 2  (nstages, ndx)     // ki, dx_rk, dli_dx
            // 3  (nstages, nx)      // y
            // 4  (nstages, nw)      // ws
            // 5  (nstages, ndx, ndx)// dki_dx, dyi_dx, ddli_ddx, Lxx_partialx
            // 6  (nstages, ndx, nu) // dki_du, dyi_du, ddli_dxdu, Lxu_i, Lxx_partialu
            // 7  (nstages, nu)      // dli_du
            // 8  (nstages, nw, nw)  // ddli_ddw
            // 9  (nstages, ndx, nw) // ddli_dxdw
            // 10 (nstages, nw, nu)  // ddli_dwdu
            // 11 (nstages, nu, nu)  // ddli_ddu, Luu_partialx

            /* ---------------------------------------------------------------- */
            /* Meta template types */
            /* ---------------------------------------------------------------- */
            using State_t = StateTpl<PhaseSpec>;

            using ActuationMeta_t = ActuationTpl<PhaseSpec>;
            using ActuationModel_t = typename ActuationMeta_t::Model;
            using ActuationData_t = typename ActuationMeta_t::Data;

            using ConstraintManagerMeta_t = ConstraintManagerTpl<PhaseSpec>;
            using ConstraintCollection_t = typename ConstraintManagerMeta_t::Collection;
            using ConstraintManagerModel_t = typename ConstraintManagerMeta_t::Model;
            using ConstraintManagerData_t = typename ConstraintManagerMeta_t::Data;

            using CostManagerMeta_t = CostManagerTpl<PhaseSpec>;
            using CostCollection_t = typename CostManagerMeta_t::Collection;
            using CostManagerModel_t = typename CostManagerMeta_t::Model;
            using CostManagerData_t = typename CostManagerMeta_t::Data;

            using NodeMeta_t = NodeTpl<PhaseSpec>;
            using NodeModel_t = typename NodeMeta_t::Model;
            using NodeData_t = typename NodeMeta_t::Data;
            using NodeDataVector_t = std::vector<NodeData_t>;

            using ControlParamMeta_t = ControlParamTpl<PhaseSpec>;
            using ControlParamModel_t = typename ControlParamMeta_t::Model;
            using ControlParamData_t = typename ControlParamMeta_t::Data;
            using ControlParamDataVector_t = std::vector<ControlParamData_t>;

            using SegmentMeta_t = SegmentTpl<PhaseSpec>;
            using SegmentModel_t = typename SegmentMeta_t::Model;
            using SegmentData_t = typename SegmentMeta_t::Data;
            using SegmentDataVector_t = std::vector<SegmentData_t>;

            using PhaseMeta_t = PhaseTpl<PhaseSpec>;
            using PhaseModel_t = typename PhaseMeta_t::Model;
            using PhaseData_t = typename PhaseMeta_t::Data;
            using PhaseDataVector_t = std::vector<PhaseData_t>;
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_spec_hpp__
