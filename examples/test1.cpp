#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "galileo/predictive/phases/phase-spec.hpp"

#include "galileo/domains/multibody/core/actuations/impl/actuation-floating-base.hpp"
#include "galileo/core/actuations/impl/actuation-full.hpp"

#include "galileo/core/states/impl/state-euclidean.hpp"
#include "galileo/domains/multibody/core/states/impl/state-multibody.hpp"

#include "galileo/core/activations/impl/activation-quadratic.hpp"
#include "galileo/core/activations/impl/activation-weighted-quadratic.hpp"

#include "galileo/core/residuals/impl/residual-control.hpp"

#include "galileo/domains/multibody/core/residuals/impl/residual-com-position.hpp"
#include "galileo/domains/multibody/core/residuals/impl/residual-frame-translation.hpp"
#include "galileo/domains/multibody/core/residuals/impl/residual-frame-velocity.hpp"
#include "galileo/domains/multibody/core/residuals/impl/residual-multibody-state.hpp"

#include "galileo/core/costs/cost-manager.hpp"
#include "galileo/core/costs/impl/cost-residual.hpp"

#include "galileo/core/constraints/equality/constraint-manager.hpp"
#include "galileo/core/constraints/equality/impl/constraint-residual.hpp"

#include "galileo/domains/multibody/spatial/contacts/impl/contact-3d.hpp"
#include "galileo/domains/multibody/spatial/contacts/impl/contact-6d.hpp"

#include "galileo/domains/multibody/spatial/contacts/contact-manager.hpp"

#include "galileo/domains/multibody/spatial/impulses/impl/impulse-3d.hpp"
#include "galileo/domains/multibody/spatial/impulses/impl/impulse-6d.hpp"

#include "galileo/domains/multibody/spatial/impulses/impulse-manager.hpp"

#include "galileo/domains/multibody/predictive/nodes/impl/node-contact-fwddyn.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include "galileo/core/controls/impl/control-param-polynomial.hpp"

#include "galileo/predictive/segments/impl/segment-erk-euler.hpp"

#include "galileo/domains/multibody/predictive/jumps/impl/jump-impulse-fwddyn.hpp"

#include "galileo/predictive/phases/impl/phase-default.hpp"

#include "galileo/core/data/data-collector-default.hpp"

#include "galileo/predictive/phases/smart-visitors/fold-engine.hpp"

#include <iostream>
#include <string>
#include <cassert>
#include <chrono>

using VarScalar = double;
using NumScalar = double;
constexpr int Options = Eigen::ColMajor;

using BasicSpec_t = galileo::BasicSpecTpl<VarScalar, NumScalar, Options>;

template <typename RobotSpec>
using StateTpl = galileo::StateMultibodyTpl<RobotSpec>;

template <typename RobotSpec>
using ActuationTpl = galileo::ActuationFloatingBaseTpl<RobotSpec>;

constexpr int NQb = 7;
constexpr int NQj = 12;
constexpr int NVb = 6;
constexpr int NVj = 12;
constexpr int NRotors = 0;

using RobotSpec_t = galileo::RobotSpecTpl<BasicSpec_t, NQb, NQj, NVb, NVj, NRotors, StateTpl, ActuationTpl>;

template <typename PhaseSpec>
using ResidualTestTpl = galileo::ResidualFrameTranslationTpl<PhaseSpec>;

template <typename PhaseSpec>
using ActivationTestTpl = galileo::ActivationQuadraticTpl<PhaseSpec, ResidualTestTpl>;

template <typename PhaseSpec>
using CostTestTpl = galileo::CostResidualTpl<PhaseSpec, ResidualTestTpl, ActivationTestTpl>;
template <typename PhaseSpec>
using CostModelTestTpl = typename galileo::traits<CostTestTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostDataTestTpl = typename galileo::traits<CostTestTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct CostCollectionTestTpl
{
    using PS = PhaseSpec;
    using CostModelVariant_t = boost::variant<CostModelTestTpl<PS>>;
    using CostDataVariant_t = boost::variant<CostDataTestTpl<PS>>;
}; // struct CostCollectionTestTpl

template <typename PhaseSpec>
using ConstraintTestTpl = galileo::ConstraintResidualTpl<PhaseSpec, ResidualTestTpl>;
template <typename PhaseSpec>
using ConstraintModelTestTpl = typename galileo::traits<ConstraintTestTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ConstraintDataTestTpl = typename galileo::traits<ConstraintTestTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ConstraintCollectionTestTpl
{
    using PS = PhaseSpec;
    using ConstraintModelVariant_t = boost::variant<ConstraintModelTestTpl<PS>>;
    using ConstraintDataVariant_t = boost::variant<ConstraintDataTestTpl<PS>>;
}; // struct ConstraintCollectionTestTpl

template <typename PhaseSpec>
using ContactTestTpl = galileo::Contact3dTpl<PhaseSpec>;
template <typename PhaseSpec>
using ContactModelTestTpl = typename galileo::traits<ContactTestTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ContactDataTestTpl = typename galileo::traits<ContactTestTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ContactCollectionTestTpl
{
    using PS = PhaseSpec;
    using ContactModelVariant_t = boost::variant<ContactModelTestTpl<PS>>;
    using ContactDataVariant_t = boost::variant<ContactDataTestTpl<PS>>;
}; // struct ContactCollectionTestTpl

template <typename PhaseSpec>
using ImpulseTestTpl = galileo::Impulse3dTpl<PhaseSpec>;
template <typename PhaseSpec>
using ImpulseModelTestTpl = typename galileo::traits<ImpulseTestTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ImpulseDataTestTpl = typename galileo::traits<ImpulseTestTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ImpulseCollectionTestTpl
{
    using PS = PhaseSpec;
    using ImpulseModelVariant_t = boost::variant<ImpulseModelTestTpl<PS>>;
    using ImpulseDataVariant_t = boost::variant<ImpulseDataTestTpl<PS>>;
}; // struct ImpulseCollectionTestTpl

template <typename PhaseSpec>
using ConstraintManagerTestTpl = galileo::ConstraintManagerTpl<PhaseSpec, ConstraintCollectionTestTpl>;

template <typename PhaseSpec>
using CostManagerTestTpl = galileo::CostManagerTpl<PhaseSpec, CostCollectionTestTpl>;

template <typename PhaseSpec>
using ContactManagerTestTpl = galileo::ContactManagerTpl<PhaseSpec, ContactCollectionTestTpl>;
template <typename PhaseSpec>
using ContactModelManagerTestTpl = galileo::ContactModelManagerTpl<PhaseSpec, ContactCollectionTestTpl>;

template <typename PhaseSpec>
using ImpulseManagerTestTpl = galileo::ImpulseManagerTpl<PhaseSpec, ImpulseCollectionTestTpl>;
template <typename PhaseSpec>
using ImpulseModelManagerTestTpl = galileo::ImpulseModelManagerTpl<PhaseSpec, ImpulseCollectionTestTpl>;

template <typename PhaseSpec>
using NodeTpl = galileo::NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTestTpl>;

static constexpr int NOrder = 1;
template <typename PhaseSpec>
using ControlParamTpl = galileo::ControlParamPolynomialTpl<PhaseSpec, NOrder>;

template <typename PhaseSpec>
using SegmentTpl = galileo::SegmentERKEulerTpl<PhaseSpec>;

template <typename PhaseSpec>
using JumpTpl = galileo::JumpImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionTestTpl>;

template <typename PhaseSpec>
using PhaseTpl = galileo::PhaseDefaultTpl<PhaseSpec, JumpTpl>;

// template <typename PhaseSpec>
// struct DummyPhaseTpl;
// template <typename PhaseSpec>
// struct DummyPhaseModelTpl;
// template <typename PhaseSpec>
// struct DummyPhaseDataTpl;

// namespace galileo
// {
//     template <typename PhaseSpec>
//     struct traits<DummyPhaseTpl<PhaseSpec>>
//     {
//         using PS = PhaseSpec;
//         using Meta_t = DummyPhaseTpl<PS>;
//         using Model_t = DummyPhaseModelTpl<PS>;
//         using Data_t = DummyPhaseDataTpl<PS>;
//     }; // traits<DummyPhaseTpl<PhaseSpec>>
// } // namespace galileo

// template <typename PhaseSpec>
// using PhaseTpl = DummyPhaseTpl<PhaseSpec>;

using PhaseSpec_t = galileo::PhaseSpecTpl<RobotSpec_t, ConstraintManagerTestTpl, CostManagerTestTpl, NodeTpl, ControlParamTpl, SegmentTpl, PhaseTpl>;

using RobotSpec_t = typename PhaseSpec_t::RS;

using RobotModel_t = typename PhaseSpec_t::RobotModel_t;
using RobotData_t = typename PhaseSpec_t::RobotData_t;

using State_t = typename PhaseSpec_t::State_t;
using ActuationModel_t = typename PhaseSpec_t::ActuationModel_t;

using ResidualModel_t = typename galileo::traits<ResidualTestTpl<PhaseSpec_t>>::Model_t;

using ActivationModel_t = typename galileo::traits<ActivationTestTpl<PhaseSpec_t>>::Model_t;

using CostModel_t = typename galileo::traits<CostTestTpl<PhaseSpec_t>>::Model_t;
using CostModelManager_t = typename PhaseSpec_t::CostModelManager_t;

using ConstraintModel_t = typename galileo::traits<ConstraintTestTpl<PhaseSpec_t>>::Model_t;
using ConstraintModelManager_t = typename PhaseSpec_t::ConstraintModelManager_t;

using ContactModel_t = typename galileo::traits<ContactTestTpl<PhaseSpec_t>>::Model_t;
using ContactModelManager_t = ContactModelManagerTestTpl<PhaseSpec_t>;

using ImpulseModel_t = typename galileo::traits<ImpulseTestTpl<PhaseSpec_t>>::Model_t;
using ImpulseModelManager_t = ImpulseModelManagerTestTpl<PhaseSpec_t>;

using NodeModel_t = typename PhaseSpec_t::NodeModel_t;
using NodeData_t = typename PhaseSpec_t::NodeData_t;

using JacobiRoots_t = galileo::JacobiRootsTpl<VarScalar, NOrder, Options>;
using BarycentricInterpolator_t = galileo::BarycentricInterpolatorTpl<VarScalar, NOrder, Options>;

using ControlParamModel_t = typename PhaseSpec_t::ControlParamModel_t;

using SegmentModel_t = typename PhaseSpec_t::SegmentModel_t;
using SegmentData_t = typename PhaseSpec_t::SegmentData_t;

using JumpModel_t = typename galileo::traits<JumpTpl<PhaseSpec_t>>::Model_t;
using JumpData_t = typename galileo::traits<JumpTpl<PhaseSpec_t>>::Data_t;

struct TestInteriorPropagator;
namespace galileo
{
    template <>
    struct traits<TestInteriorPropagator>
    {
        using FoldStateType = Eigen::VectorXd;
        using ReturnType = Eigen::VectorXd;
    };
}

// Example custom interior propagator
struct TestInteriorPropagator
    : galileo::fusion::InteriorPropagatorBase<TestInteriorPropagator>
{
    using FoldStateType = typename galileo::traits<TestInteriorPropagator>::FoldStateType;
    using ReturnType = typename galileo::traits<TestInteriorPropagator>::ReturnType;
    using ArgsType = boost::fusion::vector<const Eigen::VectorXd &>;

    // Version with control parameters
    template <typename SegmentModel, typename SegmentData>
    static ReturnType algo(const SegmentModel &segment_model, SegmentData &segment_data, FoldStateType state, const Eigen::VectorXd &controls)
    {
        std::cout << "Processing segment with controls..." << std::endl;

        // Apply segment dynamics with controls
        segment_model.calc(segment_data, state, controls);
        std::cout << "  State norm changed from " << state.norm() << " to " << segment_data.XNext_accessor().norm() << std::endl;
        segment_model.calcDiff(segment_data, state, controls);

        return segment_data.XNext_accessor();
    }
};

template <typename SegmentModel, typename SegmentData>
inline typename TestInteriorPropagator::ReturnType interior_propagate(const SegmentModel &segment_model, SegmentData &segment_data, typename TestInteriorPropagator::FoldStateType state, const Eigen::VectorXd &controls)
{
    using Algo = TestInteriorPropagator;
    return Algo::run(segment_model, segment_data, state, typename Algo::ArgsType(controls));
}

struct TestBoundaryPropagator;
namespace galileo
{
    template <>
    struct traits<TestBoundaryPropagator>
    {
        using FoldStateType = Eigen::VectorXd;
        static constexpr bool IsDirectionallyInvariant = false;
        using ReturnType = Eigen::VectorXd;
    };
}

// Example custom boundary propagator
struct TestBoundaryPropagator
    : galileo::fusion::BoundaryPropagatorBase<TestBoundaryPropagator>
{
    using FoldStateType = typename galileo::traits<TestBoundaryPropagator>::FoldStateType;
    static constexpr bool IsDirectionallyInvariant = galileo::traits<TestBoundaryPropagator>::IsDirectionallyInvariant;
    using ReturnType = typename galileo::traits<TestBoundaryPropagator>::ReturnType;
    using ArgsType = boost::fusion::vector<const Eigen::VectorXd &>;

    template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData>
    static ReturnType algo(const CurrentPhaseModel &current_phase_model, const CurrentPhaseData &current_phase_data,
                           const NextPhaseModel &next_phase_model, const NextPhaseData &next_phase_data,
                           FoldStateType state, const Eigen::VectorXd &controls)
    {
        std::cout << "Applying boundary transformation between phases" << std::endl;
        std::cout << "  Input state norm: " << state.norm() << std::endl;

        // For this test, apply a simple scaling as a "reset map"
        ReturnType reset_state = state * 0.99; // Slight energy dissipation

        std::cout << "  Output state norm: " << reset_state.norm() << std::endl;

        return reset_state;
    }
};

template <typename CurrentPhaseModel, typename CurrentPhaseData, typename NextPhaseModel, typename NextPhaseData>
inline typename TestBoundaryPropagator::ReturnType boundary_propagate(const CurrentPhaseModel &current_phase_model, const CurrentPhaseData &current_phase_data,
                                                                      const NextPhaseModel &next_phase_model, const NextPhaseData &next_phase_data,
                                                                      typename TestBoundaryPropagator::FoldStateType state, const Eigen::VectorXd &controls)
{
    using Algo = TestBoundaryPropagator;
    return Algo::run(current_phase_model, current_phase_data, next_phase_model, next_phase_data, state, typename Algo::ArgsType(controls));
}

// Test propagators for both directions
using TestPropagatorLeftFold = galileo::fusion::FoldTpl<
    TestInteriorPropagator,
    TestBoundaryPropagator,
    true>; // Left fold (forward)

using TestPropagatorRightFold = galileo::fusion::FoldTpl<
    TestInteriorPropagator,
    TestBoundaryPropagator,
    false>; // Right fold (backward)

int main()
{
    std::cout << "=== Sequential Propagation Test ===" << std::endl;

    // Setup robot model (simplified from go1_trot.cpp)
    std::string urdf_path = "/home/quant/research/Galileo/resources/go1/urdf/go1.urdf";

    RobotModel_t model = RobotModel_t();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    State_t state = State_t(model);
    ActuationModel_t actuation = ActuationModel_t(state);

    PhaseSpec_t ps = PhaseSpec_t(state);

    std::cout << "Robot initialized with " << model.nq << " positions and " << model.nv << " velocities" << std::endl;

    pinocchio::FrameIndex frame_id = model.getFrameId("base");
    ResidualModel_t residual = ResidualModel_t(ps, frame_id, Eigen::Vector3d(0., 0., 0.));
    ActivationModel_t activation = ActivationModel_t(ps, residual.get_nr_dim());

    std::cout << "Residual dimension: " << residual.get_nr_dim() << std::endl;

    CostModel_t cost = CostModel_t(ps, residual, activation);
    std::cout << "Cost model created" << std::endl;
    CostModelManager_t cost_manager = CostModelManager_t(ps);
    std::cout << "Cost model manager created" << std::endl;
    cost_manager.addItem("test_cost", cost, 1.0);
    std::cout << "Cost model added" << std::endl;

    ConstraintModel_t constraint = ConstraintModel_t(ps, residual);
    ConstraintModelManager_t constraint_manager = ConstraintModelManager_t(ps);
    constraint_manager.addItem("test_constraint", constraint);

    std::cout << "Constraint model created" << std::endl;

    ContactModelManager_t contact_manager = ContactModelManager_t(ps);
    pinocchio::FrameIndex frame_id_1 = model.getFrameId("LF_FOOT");
    ContactModel_t contact1 = ContactModel_t(ps, frame_id_1, pinocchio::LOCAL, Eigen::Vector3d(0., 0., 0.), Eigen::Vector2d(0., 0.));
    std::cout << "Contact model 1 created" << std::endl;
    contact_manager.addItem("test_contact1", contact1);
    std::cout << "Contact model 1 added" << std::endl;

    pinocchio::FrameIndex frame_id_2 = model.getFrameId("RF_FOOT");
    ContactModel_t contact2 = ContactModel_t(ps, frame_id_2, pinocchio::LOCAL, Eigen::Vector3d(0., 0., 0.), Eigen::Vector2d(0., 0.));
    std::cout << "Contact model 2 created" << std::endl;
    contact_manager.addItem("test_contact2", contact2);
    std::cout << "Contact model 2 added" << std::endl;

    // Setup basic node and control for segments
    NodeModel_t node = NodeModel_t(
        ps,
        cost_manager,
        constraint_manager,
        contact_manager,
        actuation,
        0.0,
        false);
    std::cout << "Node model created" << std::endl;

    JacobiRoots_t jacobi_roots(1.0, 0.0);
    jacobi_roots.compute_roots();
    Eigen::VectorXd nodes = jacobi_roots.get_roots();

    std::cout << "Jacobi roots created" << std::endl;

    BarycentricInterpolator_t interpolator(nodes);
    ControlParamModel_t control_param(ps, interpolator);
    const NumScalar period = 0.1;

    std::cout << "ps: " << ps << std::endl;

    // // Create test segments
    // std::vector<SegmentModel_t *> test_segments;
    // std::vector<SegmentData_t *> test_segment_data;

    // SegmentModel_t segment_model1 = SegmentModel_t(ps, node, control_param, period);
    // std::cout << "segment_model1 created" << std::endl;
    // test_segments.push_back(&segment_model1);
    // std::cout << "segment_model1 pushed" << std::endl;
    // SegmentData_t segment_data1 = segment_model1.createData();
    // std::cout << "segment_data1 created" << std::endl;
    // test_segment_data.push_back(&segment_data1);
    // std::cout << "segment_data1 pushed" << std::endl;

    // SegmentModel_t segment_model2 = SegmentModel_t(ps, node, control_param, period);
    // std::cout << "segment_model2 created" << std::endl;
    // test_segments.push_back(&segment_model2);
    // std::cout << "segment_model2 pushed" << std::endl;
    // SegmentData_t segment_data2 = segment_model2.createData();
    // std::cout << "segment_data2 created" << std::endl;
    // test_segment_data.push_back(&segment_data2);
    // std::cout << "segment_data2 pushed" << std::endl;

    // std::cout << "Created " << test_segments.size() << " test segments" << std::endl;

    // // Create initial state
    // Eigen::VectorXd initial_state = state.rand();
    // std::cout << "Initial state size: " << initial_state.size() << std::endl;
    // std::cout << "Initial state norm: " << initial_state.norm() << std::endl;

    // // Test individual segment propagation
    // std::cout << "\n--- Testing Individual Segment Transformation ---" << std::endl;

    // Eigen::VectorXd current_state = initial_state;
    // for (size_t i = 0; i < test_segments.size(); ++i)
    // {
    //     std::cout << "Segment " << i << ":" << std::endl;
    //     Eigen::VectorXd controls = Eigen::VectorXd::Random(test_segments[i]->get_ps().get_nw_dim());
    //     current_state = interior_propagate(*test_segments[i], *test_segment_data[i], current_state, controls);
    // }

    // std::cout << "Final state after individual processing: " << current_state.norm() << std::endl;

    // // Test directional fold differences
    // std::cout << "\n--- Testing Left vs Right Fold ---" << std::endl;

    // std::cout << "Left fold (forward) processing would traverse: segments[0] -> segments[" << test_segments.size() - 1 << "]" << std::endl;
    // std::cout << "Right fold (backward) processing would traverse: segments[" << test_segments.size() - 1 << "] -> segments[0]" << std::endl;

    // // Reset state for comparison
    // Eigen::VectorXd left_fold_state = initial_state;
    // Eigen::VectorXd right_fold_state = initial_state;

    // // Process with left fold (forward)
    // for (size_t i = 0; i < test_segments.size(); ++i)
    // {
    //     Eigen::VectorXd controls = Eigen::VectorXd::Random(test_segments[i]->get_ps().get_nw_dim());
    //     left_fold_state = interior_propagate(*test_segments[i], *test_segment_data[i], left_fold_state, controls);
    // }

    // // Process with right fold (backward)
    // for (size_t i = test_segments.size(); i-- > 0;)
    // {
    //     Eigen::VectorXd controls = Eigen::VectorXd::Random(test_segments[i]->get_ps().get_nw_dim());
    //     right_fold_state = interior_propagate(*test_segments[i], *test_segment_data[i], right_fold_state, controls);
    // }

    // std::cout << "Left fold final state norm: " << left_fold_state.norm() << std::endl;
    // std::cout << "Right fold final state norm: " << right_fold_state.norm() << std::endl;
    // std::cout << "Difference in norms: " << std::abs(left_fold_state.norm() - right_fold_state.norm()) << std::endl;

    // pinocchio::FrameIndex frame_id_3 = model.getFrameId("LF_FOOT");
    // ImpulseModel_t impulse = ImpulseModel_t(ps, frame_id_3, pinocchio::LOCAL);
    // ImpulseModelManager_t impulse_manager = ImpulseModelManager_t(ps);
    // impulse_manager.addItem("test_impulse", impulse);
    // std::cout << "Impulse model created" << std::endl;

    // JumpModel_t jump = JumpModel_t(ps, cost_manager, constraint_manager, impulse_manager);
    // std::cout << "Jump model created" << std::endl;
    // JumpData_t jump_data = jump.createData();
    // std::cout << "Jump data created" << std::endl;
    // jump.calc(jump_data, initial_state);
    // std::cout << "Jump calc done" << std::endl;
    // jump.calcDiff(jump_data, initial_state);
    // std::cout << "Jump calcDiff done" << std::endl;

    // --------------------------------------------------------------------------------------------------------

    std::vector<SegmentModel_t> segment_models;
    std::vector<SegmentData_t> segment_datas;
    std::vector<JumpModel_t> jump_models;
    std::vector<JumpData_t> jump_datas;

    pinocchio::FrameIndex frame_id_3 = model.getFrameId("LF_FOOT");
    ImpulseModel_t impulse = ImpulseModel_t(ps, frame_id_3, pinocchio::LOCAL);
    ImpulseModelManager_t impulse_manager = ImpulseModelManager_t(ps);
    impulse_manager.addItem("test_impulse", impulse);
    std::cout << "Impulse model created" << std::endl;

    for (int i = 0; i < 10; ++i)
    {
        segment_models.push_back(SegmentModel_t(ps, node, control_param, period));
        segment_datas.push_back(segment_models.back().createData());
        jump_models.push_back(JumpModel_t(ps, cost_manager, constraint_manager, ImpulseModelManager_t(ps)));
        jump_datas.push_back(jump_models.back().createData());
    }

    using VectorNx_t = typename PhaseSpec_t::VectorNx_t;
    using VectorNw_t = typename PhaseSpec_t::VectorNw_t;

    for (int i = 0; i < 10; ++i)
    {
        VectorNx_t x = state.rand();
        VectorNw_t w = VectorNw_t::Random(ps.get_nw_dim());
        std::cout << "Calling segment model calc" << std::endl;
        segment_models[i].calc(segment_datas[i], x, w);
        std::cout << "Calling segment model calcDiff" << std::endl;
        segment_models[i].calcDiff(segment_datas[i], x, w);

        std::cout << "Calling jump model calc" << std::endl;
        jump_models[i].calc(jump_datas[i], x);
        std::cout << "Calling jump model calcDiff" << std::endl;
        jump_models[i].calcDiff(jump_datas[i], x);
    }

    return 0;
}
