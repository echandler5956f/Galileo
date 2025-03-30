#ifndef __galileo_predictive_phases_phase_spec_hpp__
#define __galileo_predictive_phases_phase_spec_hpp__

#include "galileo/predictive/phases/fwd.hpp"
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
            using PhaseSpec = PhaseSpecTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX, _NStages, _NOrder, StateTpl, ActuationTpl, ConstraintManagerTpl, CostManagerTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl, _SegmentType>;

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
            static constexpr int NStages = _NStages; // Number of Runge-Kutta stages per segment
            static constexpr int NOrder = _NOrder;   // Order of the control parameterization per segment
            static constexpr int NW = NU * NOrder;   // Number of control parameters

            static constexpr SegmentType SegmentType = _SegmentType;

            /* ---------------------------------------------------------------- */
            /* Node type definitions */
            /* ---------------------------------------------------------------- */
            // Dynamics
            using F_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;   // Dynamics vector
            using Fx_t = Eigen::Matrix<VarScalar, NDX, NX, Options>; // Jacobian of dynamics w.r.t. state
            using Fu_t = Eigen::Matrix<VarScalar, NDX, NU, Options>; // Jacobian of dynamics w.r.t. control

            // Cost (CostManager holds an Eigen map to these, which are stored in NodeData)
            using L_t = Eigen::Matrix<VarScalar, 1, 1, Options>;     // Cost scalar
            using Lx_t = Eigen::Matrix<VarScalar, 1, NX, Options>;   // Jacobian of cost w.r.t. state
            using Lu_t = Eigen::Matrix<VarScalar, 1, NU, Options>;   // Jacobian of cost w.r.t. control
            using Lxx_t = Eigen::Matrix<VarScalar, NX, NX, Options>; // Hessian of cost w.r.t. state
            using Lxu_t = Eigen::Matrix<VarScalar, NX, NU, Options>; // Hessian of cost w.r.t. state and control
            using Luu_t = Eigen::Matrix<VarScalar, NU, NU, Options>; // Hessian of cost w.r.t. control

            // Equality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using H_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;   // Equality constraint vector
            using Hx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NX, Options>; // Jacobian of equality constraints w.r.t. state
            using Hu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>; // Jacobian of equality constraints w.r.t. control

            // Inequality constraints (ConstraintManager holds an Eigen map to these, which are stored in NodeData)
            using G_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1, Options>;   // Inequality constraint vector
            using Gx_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NX, Options>; // Jacobian of inequality constraints w.r.t. state
            using Gu_t = Eigen::Matrix<VarScalar, Eigen::Dynamic, NU, Options>; // Jacobian of inequality constraints w.r.t. control

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

            /* ---------------------------------------------------------------- */
            /* An assortment of Eigen types (primarily for use in Segments) */
            /* ---------------------------------------------------------------- */
            using VectorNdx = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using VectorNdxArray = std::array<Eigen::Matrix<VarScalar, NDX, 1, Options>, NStages>;
            using VectorNxArray = std::array<Eigen::Matrix<VarScalar, NX, 1, Options>, NStages>;
            using VectorNuArray = std::array<Eigen::Matrix<VarScalar, NU, 1, Options>, NStages>;
            using MatrixNdxArray = std::array<Eigen::Matrix<VarScalar, NDX, NDX, Options>, NStages>;
            using MatrixNdxNwArray = std::array<Eigen::Matrix<VarScalar, NDX, NW, Options>, NStages>;
            using VectorNwArray = std::array<Eigen::Matrix<VarScalar, NW, 1, Options>, NStages>;
            using MatrixNuArray = std::array<Eigen::Matrix<VarScalar, NU, NU, Options>, NStages>;
            using MatrixNdxNuArray = std::array<Eigen::Matrix<VarScalar, NDX, NU, Options>, NStages>;
            using MatrixNuNwArray = std::array<Eigen::Matrix<VarScalar, NU, NW, Options>, NStages>;
            using MatrixNwArray = std::array<Eigen::Matrix<VarScalar, NW, NW, Options>, NStages>;

            /* ---------------------------------------------------------------- */
            /* Meta template types */
            /* ---------------------------------------------------------------- */
            using State_t = StateTpl<PhaseSpec>;

            using ActuationMeta_t = ActuationTpl<PhaseSpec>;
            using ActuationModel_t = typename ActuationMeta_t::Model;
            using ActuationData_t = typename ActuationMeta_t::Data;

            using ConstraintManagerMeta_t = ConstraintManagerTpl<PhaseSpec>;
            using ConstraintManagerModel_t = typename ConstraintManagerMeta_t::Model;
            using ConstraintManagerData_t = typename ConstraintManagerMeta_t::Data;

            using CostManagerMeta_t = CostManagerTpl<PhaseSpec>;
            using CostManagerModel_t = typename CostManagerMeta_t::Model;
            using CostManagerData_t = typename CostManagerMeta_t::Data;

            using NodeMeta_t = NodeTpl<PhaseSpec>;
            using NodeModel_t = typename NodeMeta_t::Model;
            using NodeData_t = typename NodeMeta_t::Data;

            using ControlParamMeta_t = ControlParamTpl<PhaseSpec>;
            using ControlParamModel_t = typename ControlParamMeta_t::Model;
            using ControlParamData_t = typename ControlParamMeta_t::Data;

            using SegmentMeta_t = SegmentTpl<PhaseSpec>;
            using SegmentModel_t = typename SegmentMeta_t::Model;
            using SegmentData_t = typename SegmentMeta_t::Data;

            using PhaseMeta_t = PhaseTpl<PhaseSpec>;
            using PhaseModel_t = typename PhaseMeta_t::Model;
            using PhaseData_t = typename PhaseMeta_t::Data;
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_spec_hpp__
