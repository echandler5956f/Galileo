// #include <galileo/galileo.h>

#include <galileo/fwd.hpp>
#include <galileo/core/states/multibody.hpp>
#include <galileo/core/actuations/floating-base.hpp>
#include <galileo/core/basic-spec.hpp>
#include <galileo/predictive/phases/phase-spec.hpp>

#include <pinocchio/multibody/model.hpp>

using namespace galileo::core;
using namespace galileo::predictive;

// These template classes basically create a code generator for the optimal control problem by means of mixins
using BasicSpec = BasicSpecTpl<double, double, Eigen::RowMajor, 7, 12, 6, 12, 0, StateMultibodyTpl, ActuationFloatingBaseTpl>;

using ConstraintManagerTpl = ConstraintManagerDefaultTpl;
using CostManagerTpl = CostManagerDefaultTpl;
using ContactManagerTpl = ContactManagerDefaultTpl;

using NodeTpl = NodeContactFwdDynTpl<ContactManagerTpl>;

static constexpr int NOrder = 2;
using ControlParamTpl = ControlParamJacobiPolynomialTpl<NOrder>;

using SegmentTpl = SegmentERKEulerTpl;

using ImpulseManagerTpl = ImpulseManagerDefaultTpl;

using PhaseTpl = PhaseWithImpulsePolicyTpl<ImpulseManagerTpl>;

using PhaseSpec = PhaseSpecTpl<BasicSpec, ConstraintManagerTpl, CostManagerTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl>;

using Phase1 = PhaseTpl<PhaseSpec>;
using Phase2 = PhaseTpl<PhaseSpec>;

using PhaseCollection_t = PhaseCollectionTpl<Phase1, Phase2>;

using Trajectory_t = TrajectoryTpl<PhaseCollection_t>;

using OCP_t = OCPTpl<Trajectory_t>;

using Solver_t = SolverDDPTpl<OCP_t>;

int main(int argc, char *argv[])
{
    // Create a quadruped walking problem
    // 1. Define the phase specs
    // 2. Define the phase collection
    // 3. Define the trajectory
    // 4. Define the OCP
    // 5. Solve the OCP

    return 0;
}
