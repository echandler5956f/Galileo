#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <galileo/fwd.hpp>
#include <galileo/core/robot-spec.hpp>

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

// #include <galileo/predictive/trajectory.hpp>

// #include <galileo/predictive/optimal-control-problem.hpp>

// #include <galileo/predictive/solvers/ddp.hpp>

#include <iostream>

using namespace galileo;

template <typename RS>
using StateTpl = StateMultibodyTpl<RS>;

template <typename RS>
using ActuationTpl = ActuationFloatingBaseTpl<RS>;

constexpr int NQb = 7;
constexpr int NQj = 12;
constexpr int NVb = 6;
constexpr int NVj = 12;
constexpr int NRotors = 0;

using VarScalar = double;
using NumScalar = double;
constexpr int Options = Eigen::ColMajor;

// These template classes basically create a code generator for the optimal control problem by means of mixins
using RobotSpec_t = RobotSpecTpl<VarScalar, NumScalar, Options, NQb, NQj, NVb, NVj, NRotors, StateTpl, ActuationTpl>;

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
template <typename PS_>
struct DummyPhaseModelTpl;
template <typename PS_>
struct DummyPhaseDataTpl;

using PhaseSpec_t = PhaseSpecTpl<RobotSpec_t, ConstraintManagerDefaultTpl, CostManagerDefaultTpl, NodeTpl, ControlParamTpl, SegmentTpl, DummyPhaseTpl>;

// template <typename tmpScalar1, typename tmpScalar2, int tmpOptions1>
// using PhaseSpecDummyTpl = PhaseSpec_t;

// template <typename PhaseSpec1, typename PhaseSpec2>
// struct PhaseCollectionExample
// {
// public:

//     using Phase1Model_t = DummyPhaseModelTpl<PhaseSpec1>;
//     using Phase2Model_t = DummyPhaseModelTpl<PhaseSpec2>;

//     using Phase1Data_t = DummyPhaseDataTpl<PhaseSpec1>;
//     using Phase2Data_t = DummyPhaseDataTpl<PhaseSpec2>;

//     using ModelVariant_t = boost::variant<Phase1Model_t, Phase2Model_t>;
//     using DataVariant_t = boost::variant<Phase1Data_t, Phase2Data_t>;
// };

// template <typename tmpScalar1, typename tmpScalar2, int tmpOptions1>
// using PhaseCollectionDummyTpl = PhaseCollectionExample<PhaseSpecDummyTpl<tmpScalar1, tmpScalar2, tmpOptions1>, PhaseSpecDummyTpl<tmpScalar1, tmpScalar2, tmpOptions1>>;

// using Trajectory_t = Trajectory<VarScalar, NumScalar, Options, PhaseCollectionDummyTpl>;

// using OCP_t = OptimalControlProblem<VarScalar, NumScalar, Options, PhaseCollectionDummyTpl>;

// using Solver_t = SolverDDP<VarScalar, NumScalar, Options, FeasibilityNormOptions::L1, PhaseCollectionDummyTpl>;

int main(int argc, char *argv[])
{
    std::string urdf_path = "/home/quant/Galileo/resources/go1/urdf/go1.urdf";
    pinocchio::ModelTpl<VarScalar, Options> model = pinocchio::ModelTpl<VarScalar, Options> ();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    StateTpl<RobotSpec_t> state(&model);

    // Test the state
    std::cout << "state.zero() = " << state.zero() << std::endl;
    std::cout << "state.rand() = " << state.rand() << std::endl;

    // Test the actuation
    using ActuationModel_t = RobotSpec_t::ActuationModel_t;
    using ActuationData_t = RobotSpec_t::ActuationData_t;
    ActuationModel_t actuation;
    ActuationData_t actuation_data = actuation.createData();
    std::cout << "actuation_data.dtau_du = " << actuation_data.dtau_du << std::endl;
    std::cout << "actuation_data.Mtau = " << actuation_data.Mtau << std::endl;
    std::cout << "actuation_data.u = " << actuation_data.u << std::endl;
    std::cout << "actuation_data.tau = " << actuation_data.tau << std::endl;

    // Test the constraint manager
    
    // Create a quadruped walking problem
    // 1. Define the phase specs
    // 2. Define the phase collection
    // 3. Define the trajectory
    // 4. Define the OCP
    // 5. Solve the OCP

    return 0;
}
