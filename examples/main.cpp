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
using BasicSpec_t = BasicSpecTpl<VarScalar, NumScalar, Options, NQb, NQj, NVb, NVj, NRotors, StateTpl, ActuationTpl>;

template <typename PS_>
using ConstraintManagerDefaultTpl = ConstraintManagerTpl<PS_, ConstraintCollectionDefaultTpl>;

template <typename PS_>
using CostManagerDefaultTpl = CostManagerTpl<PS_, CostCollectionDefaultTpl>;

template <typename PS_>
using ContactManagerDefaultTpl = ContactManagerTpl<PS_, ContactCollectionDefaultTpl>;

template <typename PS_>
using NodeTpl = NodeContactFwdDynTpl<PS_, ContactManagerDefaultTpl>;

static constexpr int NOrder = 2;
template <typename PS_>
using ControlParamTpl = ControlParamJacobiPolynomialTpl<PS_, NOrder>;

template <typename PS_>
using SegmentTpl = SegmentERKEulerTpl<PS_>;

// template <typename PS_>
// using ImpulseManagerTpl = ImpulseManagerTpl<PS_>;

// template <typename PS_>
// using PhaseTpl = PhaseWithImpulsePolicyTpl<PS_, ImpulseManagerTpl>;

template <typename PS_>
struct DummyPhaseTpl;

using PhaseSpec_t = PhaseSpecTpl<BasicSpec_t, ConstraintManagerDefaultTpl, CostManagerDefaultTpl, NodeTpl, ControlParamTpl, SegmentTpl, DummyPhaseTpl>;

using Phase1_t = DummyPhaseTpl<PhaseSpec_t>;
using Phase2_t = DummyPhaseTpl<PhaseSpec_t>;

using PhaseCollection_t = PhaseCollectionTpl<Phase1_t, Phase2_t>;

using Trajectory_t = TrajectoryTpl<PhaseCollection_t>;

using OCP_t = OCPTpl<Trajectory_t>;

using Solver_t = SolverDDPTpl<OCP_t>;

int main(int argc, char *argv[])
{
    std::string urdf_path = "/home/quant/Galileo/resources/go1/urdf/go1.urdf";
    pinocchio::ModelTpl<VarScalar, Options> model = pinocchio::ModelTpl<VarScalar, Options> ();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    StateTpl<BasicSpec_t> state(&model);

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
