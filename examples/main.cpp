#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <galileo/fwd.hpp>
#include <galileo/core/basic-spec.hpp>

#include <galileo/core/states/state-base.hpp>
#include <galileo/core/states/multibody.hpp>

#include <galileo/core/actuations/actuation-base.hpp>
#include <galileo/core/actuations/floating-base.hpp>

#include <galileo/predictive/phases/phase-spec.hpp>

#include <galileo/core/constraints/constraint-manager.hpp>
#include <galileo/core/costs/cost-manager.hpp>
#include <galileo/multibody/contacts/contact-manager.hpp>

#include <galileo/predictive/nodes/node-contact-fwddyn.hpp>

#include <galileo/core/controls/control-param-jpoly.hpp>

#include <galileo/predictive/segments/segment-erk-euler.hpp>

#include <iostream>

using namespace galileo;

template <typename BS>
using StateTpl = StateMultibodyTpl<BS>;

template <typename BS>
using ActuationTpl = ActuationFloatingBaseTpl<BS>;

constexpr int NQb = 7;
constexpr int NQj = 12;
constexpr int NVb = 6;
constexpr int NVj = 12;
constexpr int NRotors = 0;

using VarScalar = double;
using NumScalar = double;
constexpr int Options = Eigen::ColMajor;

// These template classes basically create a code generator for the optimal control problem by means of mixins
using BasicSpec = BasicSpecTpl<VarScalar, NumScalar, Options, NQb, NQj, NVb, NVj, NRotors, StateTpl, ActuationTpl>;

using ConstraintManagerTpl = ConstraintManagerDefaultTpl;
using CostManagerTpl = CostManagerDefaultTpl;
using ContactManagerTpl = ContactManagerDefaultTpl;

template <typename PS>
using NodeTpl = NodeContactFwdDynTpl<PS, ContactManagerTpl>;

static constexpr int NOrder = 2;
template <typename PS>
using ControlParamTpl = ControlParamJacobiPolynomialTpl<PS, NOrder>;

template <typename PS>
using SegmentTpl = SegmentERKEulerTpl<PS>;

using ImpulseManagerTpl = ImpulseManagerDefaultTpl;

template <typename PS>
using PhaseTpl = PhaseWithImpulsePolicyTpl<PS, ImpulseManagerTpl>;

using PhaseSpec = PhaseSpecTpl<BasicSpec, ConstraintManagerTpl, CostManagerTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl>;

using Phase1 = PhaseTpl<PhaseSpec>;
using Phase2 = PhaseTpl<PhaseSpec>;

using PhaseCollection_t = PhaseCollectionTpl<Phase1, Phase2>;

using Trajectory_t = TrajectoryTpl<PhaseCollection_t>;

using OCP_t = OCPTpl<Trajectory_t>;

using Solver_t = SolverDDPTpl<OCP_t>;

int main(int argc, char *argv[])
{
    std::string urdf_path = "/home/quant/Galileo/resources/go1/urdf/go1.urdf";
    pinocchio::ModelTpl<VarScalar, Options> model = pinocchio::ModelTpl<VarScalar, Options> ();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    StateTpl<BasicSpec> state(&model);

    // Test the state
    std::cout << "state.zero() = " << state.zero() << std::endl;
    std::cout << "state.rand() = " << state.rand() << std::endl;

    // Create a quadruped walking problem
    // 1. Define the phase specs
    // 2. Define the phase collection
    // 3. Define the trajectory
    // 4. Define the OCP
    // 5. Solve the OCP

    return 0;
}
