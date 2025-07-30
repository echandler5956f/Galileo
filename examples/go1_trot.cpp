#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "galileo/predictive/phases/phase-spec.hpp"

#include "galileo/core/actuations/implementations/actuation-floating-base.hpp"

#include "galileo/core/states/implementations/state-multibody.hpp"

#include "galileo/core/activations/implementations/activation-weighted-quadratic.hpp"

#include "galileo/core/residuals/implementations/residual-com-position.hpp"
#include "galileo/core/residuals/implementations/residual-control.hpp"
#include "galileo/core/residuals/implementations/residual-frame-translation.hpp"
#include "galileo/core/residuals/implementations/residual-frame-velocity.hpp"
#include "galileo/core/residuals/implementations/residual-state.hpp"

#include "galileo/core/costs/cost-manager.hpp"
#include "galileo/core/costs/implementations/cost-residual.hpp"

#include "galileo/core/constraints/equality/constraint-manager.hpp"
#include "galileo/core/constraints/equality/implementations/constraint-residual.hpp"

#include "galileo/multibody/contacts/implementations/contact-3d.hpp"
#include "galileo/multibody/contacts/contact-manager.hpp"

#include "galileo/multibody/impulses/implementations/impulse-3d.hpp"
#include "galileo/multibody/impulses/impulse-manager.hpp"

#include "galileo/predictive/nodes/implementations/node-contact-fwddyn.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include "galileo/core/controls/implementations/control-param-polynomial.hpp"

#include "galileo/predictive/segments/implementations/segment-erk-euler.hpp"

#include "galileo/predictive/jumps/implementations/jump-impulse-fwddyn.hpp"

#include "galileo/predictive/phases/implementations/phase-default.hpp"

#include "galileo/core/data/data-collector-default.hpp"

#include "galileo/predictive/phases/phase-generic.hpp"

#include "galileo/predictive/phases/smart-visitors/fold-engine.hpp"

#include <iostream>
#include <string>

// Standard library
#include <cassert>
#include <chrono>

using VarScalar_ = double;
using NumScalar_ = double;
constexpr int Options_ = Eigen::ColMajor;

using BasicSpec_t = galileo::BasicSpecTpl<VarScalar_, NumScalar_, Options_>;

template <typename RobotSpec>
using StateTpl = galileo::StateMultibodyTpl<RobotSpec>;

template <typename RobotSpec>
using ActuationTpl = galileo::ActuationFloatingBaseTpl<RobotSpec>;

constexpr int NQb_ = 7;
constexpr int NQj_ = 12;
constexpr int NVb_ = 6;
constexpr int NVj_ = 12;
constexpr int NRotors_ = 0;

using RobotSpec_t = galileo::RobotSpecTpl<BasicSpec_t, NQb_, NQj_, NVb_, NVj_, NRotors_, StateTpl, ActuationTpl>;

// Define the 5 different residual types for walking gait
template <typename PhaseSpec>
using ResidualCoMPositionTpl = galileo::ResidualCoMPositionTpl<PhaseSpec>;
template <typename PhaseSpec>
using ResidualControlTpl = galileo::ResidualControlTpl<PhaseSpec>;
template <typename PhaseSpec>
using ResidualFrameTranslationTpl = galileo::ResidualFrameTranslationTpl<PhaseSpec>;
template <typename PhaseSpec>
using ResidualFrameVelocityTpl = galileo::ResidualFrameVelocityTpl<PhaseSpec>;
template <typename PhaseSpec>
using ResidualStateTpl = galileo::ResidualStateTpl<PhaseSpec>;

// Define activation types for each residual
template <typename PhaseSpec>
using ActivationCoMTpl = galileo::ActivationWeightedQuadraticTpl<PhaseSpec, ResidualCoMPositionTpl>;
template <typename PhaseSpec>
using ActivationControlTpl = galileo::ActivationWeightedQuadraticTpl<PhaseSpec, ResidualControlTpl>;
template <typename PhaseSpec>
using ActivationFrameTranslationTpl = galileo::ActivationWeightedQuadraticTpl<PhaseSpec, ResidualFrameTranslationTpl>;
template <typename PhaseSpec>
using ActivationFrameVelocityTpl = galileo::ActivationWeightedQuadraticTpl<PhaseSpec, ResidualFrameVelocityTpl>;
template <typename PhaseSpec>
using ActivationStateTpl = galileo::ActivationWeightedQuadraticTpl<PhaseSpec, ResidualStateTpl>;

// Define cost types for each residual
template <typename PhaseSpec>
using CostCoMPositionTpl = galileo::CostResidualTpl<PhaseSpec, ResidualCoMPositionTpl, ActivationCoMTpl>;
template <typename PhaseSpec>
using CostControlTpl = galileo::CostResidualTpl<PhaseSpec, ResidualControlTpl, ActivationControlTpl>;
template <typename PhaseSpec>
using CostFrameTranslationTpl = galileo::CostResidualTpl<PhaseSpec, ResidualFrameTranslationTpl, ActivationFrameTranslationTpl>;
template <typename PhaseSpec>
using CostFrameVelocityTpl = galileo::CostResidualTpl<PhaseSpec, ResidualFrameVelocityTpl, ActivationFrameVelocityTpl>;
template <typename PhaseSpec>
using CostStateTpl = galileo::CostResidualTpl<PhaseSpec, ResidualStateTpl, ActivationStateTpl>;

// Extract model and data types for costs
template <typename PhaseSpec>
using CostCoMPositionModelTpl = typename galileo::traits<CostCoMPositionTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostCoMPositionDataTpl = typename galileo::traits<CostCoMPositionTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
using CostControlModelTpl = typename galileo::traits<CostControlTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostControlDataTpl = typename galileo::traits<CostControlTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
using CostFrameTranslationModelTpl = typename galileo::traits<CostFrameTranslationTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostFrameTranslationDataTpl = typename galileo::traits<CostFrameTranslationTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
using CostFrameVelocityModelTpl = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostFrameVelocityDataTpl = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
using CostStateModelTpl = typename galileo::traits<CostStateTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostStateDataTpl = typename galileo::traits<CostStateTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct CostCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using CostModelVariant_t = boost::variant<CostCoMPositionModelTpl<PS>,
                                              CostControlModelTpl<PS>,
                                              CostFrameTranslationModelTpl<PS>,
                                              CostFrameVelocityModelTpl<PS>,
                                              CostStateModelTpl<PS>>;

    using CostDataVariant_t = boost::variant<CostCoMPositionDataTpl<PS>,
                                             CostControlDataTpl<PS>,
                                             CostFrameTranslationDataTpl<PS>,
                                             CostFrameVelocityDataTpl<PS>,
                                             CostStateDataTpl<PS>>;
}; // struct CostCollectionWalkingTpl

template <typename PhaseSpec>
using ConstraintTestTpl = galileo::ConstraintResidualTpl<PhaseSpec, ResidualFrameTranslationTpl>;
template <typename PhaseSpec>
using ConstraintModelTestTpl = typename galileo::traits<ConstraintTestTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ConstraintDataTestTpl = typename galileo::traits<ConstraintTestTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ConstraintCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using ConstraintModelVariant_t = boost::variant<ConstraintModelTestTpl<PS>>;
    using ConstraintDataVariant_t = boost::variant<ConstraintDataTestTpl<PS>>;
}; // struct ConstraintCollectionWalkingTpl

template <typename PhaseSpec>
using ContactWalkingTpl = galileo::Contact3dTpl<PhaseSpec>;
template <typename PhaseSpec>
using ContactModelWalkingTpl = typename galileo::traits<ContactWalkingTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ContactDataWalkingTpl = typename galileo::traits<ContactWalkingTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ContactCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using ContactModelVariant_t = boost::variant<ContactModelWalkingTpl<PS>>;
    using ContactDataVariant_t = boost::variant<ContactDataWalkingTpl<PS>>;
}; // struct ContactCollectionWalkingTpl

template <typename PhaseSpec>
using ImpulseWalkingTpl = galileo::Impulse3dTpl<PhaseSpec>;
template <typename PhaseSpec>
using ImpulseModelWalkingTpl = typename galileo::traits<ImpulseWalkingTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ImpulseDataWalkingTpl = typename galileo::traits<ImpulseWalkingTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ImpulseCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using ImpulseModelVariant_t = boost::variant<ImpulseModelWalkingTpl<PS>>;
    using ImpulseDataVariant_t = boost::variant<ImpulseDataWalkingTpl<PS>>;
}; // struct ImpulseCollectionWalkingTpl

template <typename PhaseSpec>
using ConstraintManagerWalkingTpl = galileo::ConstraintManagerTpl<PhaseSpec, ConstraintCollectionWalkingTpl>;

template <typename PhaseSpec>
using CostManagerWalkingTpl = galileo::CostManagerTpl<PhaseSpec, CostCollectionWalkingTpl>;

template <typename PhaseSpec>
using ContactManagerWalkingTpl = galileo::ContactManagerTpl<PhaseSpec, ContactCollectionWalkingTpl>;
template <typename PhaseSpec>
using ContactModelManagerWalkingTpl = galileo::ContactModelManagerTpl<PhaseSpec, ContactCollectionWalkingTpl>;

template <typename PhaseSpec>
using ImpulseManagerWalkingTpl = galileo::ImpulseManagerTpl<PhaseSpec, ImpulseCollectionWalkingTpl>;
template <typename PhaseSpec>
using ImpulseModelManagerWalkingTpl = galileo::ImpulseModelManagerTpl<PhaseSpec, ImpulseCollectionWalkingTpl>;

template <typename PhaseSpec>
using NodeWalkingTpl = galileo::NodeContactFwdDynTpl<PhaseSpec, ContactCollectionWalkingTpl>;

static constexpr int NOrder_ = 1;
template <typename PhaseSpec>
using ControlParamWalkingTpl = galileo::ControlParamPolynomialTpl<PhaseSpec, NOrder_>;

template <typename PhaseSpec>
using SegmentWalkingTpl = galileo::SegmentERKEulerTpl<PhaseSpec>;

template <typename PhaseSpec>
using JumpWalkingTpl = galileo::JumpImpulseFwdDynTpl<PhaseSpec, ImpulseCollectionWalkingTpl>;

template <typename PhaseSpec>
using PhaseWalkingTpl = galileo::PhaseDefaultTpl<PhaseSpec, JumpWalkingTpl>;

using PhaseSpec_t = galileo::PhaseSpecTpl<RobotSpec_t, ConstraintManagerWalkingTpl, CostManagerWalkingTpl, NodeWalkingTpl, ControlParamWalkingTpl, SegmentWalkingTpl, PhaseWalkingTpl>;

using RobotModel_t = typename PhaseSpec_t::RobotModel_t;
using RobotData_t = typename PhaseSpec_t::RobotData_t;

using State_t = typename PhaseSpec_t::State_t;
using ActuationModel_t = typename PhaseSpec_t::ActuationModel_t;

using CostModelManager_t = typename PhaseSpec_t::CostModelManager_t;
using ConstraintModelManager_t = typename PhaseSpec_t::ConstraintModelManager_t;
using ContactModelManager_t = ContactModelManagerWalkingTpl<PhaseSpec_t>;
using ImpulseModelManager_t = ImpulseModelManagerWalkingTpl<PhaseSpec_t>;

using NodeModel_t = typename PhaseSpec_t::NodeModel_t;
using NodeData_t = typename PhaseSpec_t::NodeData_t;

using JacobiRoots_t = galileo::JacobiRootsTpl<VarScalar_, NOrder_, Options_>;
using BarycentricInterpolator_t = galileo::BarycentricInterpolatorTpl<VarScalar_, NOrder_, Options_>;

using ControlParamModel_t = typename PhaseSpec_t::ControlParamModel_t;

using SegmentModel_t = typename PhaseSpec_t::SegmentModel_t;
using SegmentData_t = typename PhaseSpec_t::SegmentData_t;

using JumpModel_t = typename galileo::traits<JumpWalkingTpl<PhaseSpec_t>>::Model_t;
using PhaseModel_t = typename PhaseSpec_t::PhaseModel_t;

// Extract concrete model types for the walking gait
using ResidualCoMPositionModel_t = typename galileo::traits<ResidualCoMPositionTpl<PhaseSpec_t>>::Model_t;
using ActivationCoMModel_t = typename galileo::traits<ActivationCoMTpl<PhaseSpec_t>>::Model_t;
using CostCoMPositionModel_t = typename galileo::traits<CostCoMPositionTpl<PhaseSpec_t>>::Model_t;

using ResidualControlModel_t = typename galileo::traits<ResidualControlTpl<PhaseSpec_t>>::Model_t;
using ActivationControlModel_t = typename galileo::traits<ActivationControlTpl<PhaseSpec_t>>::Model_t;
using CostControlModel_t = typename galileo::traits<CostControlTpl<PhaseSpec_t>>::Model_t;

using ResidualStateModel_t = typename galileo::traits<ResidualStateTpl<PhaseSpec_t>>::Model_t;
using ActivationStateModel_t = typename galileo::traits<ActivationStateTpl<PhaseSpec_t>>::Model_t;
using CostStateModel_t = typename galileo::traits<CostStateTpl<PhaseSpec_t>>::Model_t;

using ResidualFrameTranslationModel_t = typename galileo::traits<ResidualFrameTranslationTpl<PhaseSpec_t>>::Model_t;
using ActivationFrameTranslationModel_t = typename galileo::traits<ActivationFrameTranslationTpl<PhaseSpec_t>>::Model_t;
using CostFrameTranslationModel_t = typename galileo::traits<CostFrameTranslationTpl<PhaseSpec_t>>::Model_t;

using ResidualFrameVelocityModel_t = typename galileo::traits<ResidualFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using ActivationFrameVelocityModel_t = typename galileo::traits<ActivationFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using CostFrameVelocityModel_t = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec_t>>::Model_t;

using ContactModel_t = typename galileo::traits<ContactWalkingTpl<PhaseSpec_t>>::Model_t;

template <typename BasicSpec>
struct PhaseCollectionWalkingTpl
{
    using PhaseModelVariant_t = boost::variant<PhaseModel_t>;
    using PhaseDataVariant_t = boost::variant<typename PhaseModel_t::Data_t>;
}; // struct PhaseCollectionWalkingTpl

GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PhaseSpec_t);

using Phase_t = galileo::PhaseTpl<BasicSpec_t, PhaseCollectionWalkingTpl>;
using PhaseModelGeneric_t = typename galileo::traits<Phase_t>::Model_t;
using PhaseDataGeneric_t = typename galileo::traits<Phase_t>::Data_t;

// Function to create cost manager for swing foot phase
CostModelManager_t createSwingFootCostManager(const PhaseSpec_t &ps,
                                              const Eigen::Vector3d &com_ref,
                                              const std::vector<pinocchio::FrameIndex> &swing_foot_ids = {},
                                              const std::vector<Eigen::Vector3d> &swing_foot_targets = {})
{
    CostModelManager_t cost_manager(ps);

    // CoM position cost (always active)
    ResidualCoMPositionModel_t com_residual(ps, com_ref);
    Eigen::Vector<VarScalar_, ResidualCoMPositionModel_t::DimNR_t::Value> com_activation_weight =
        Eigen::Vector<VarScalar_, ResidualCoMPositionModel_t::DimNR_t::Value>::Constant(1.0);
    ActivationCoMModel_t com_activation(ps, com_residual.get_nr_dim(), com_activation_weight);
    CostCoMPositionModel_t com_cost(ps, com_residual, com_activation);
    cost_manager.addItem("com_position", com_cost, 1e3);

    // Control regularization cost
    Eigen::Vector<VarScalar_, PhaseSpec_t::DimNU_t::Value> control_ref =
        Eigen::Vector<VarScalar_, PhaseSpec_t::DimNU_t::Value>::Zero();
    ResidualControlModel_t control_residual(ps, control_ref);
    Eigen::Vector<VarScalar_, ResidualControlModel_t::DimNR_t::Value> control_activation_weight =
        Eigen::Vector<VarScalar_, ResidualControlModel_t::DimNR_t::Value>::Constant(1.0);
    ActivationControlModel_t control_activation(ps, control_residual.get_nr_dim(), control_activation_weight);
    CostControlModel_t control_cost(ps, control_residual, control_activation);
    cost_manager.addItem("control_regularization", control_cost, 1e-1);

    // State regularization cost
    ResidualStateModel_t state_residual(ps, ps.get_state().zero());
    Eigen::Vector<VarScalar_, ResidualStateModel_t::DimNR_t::Value> state_activation_weight =
        Eigen::Vector<VarScalar_, ResidualStateModel_t::DimNR_t::Value>::Constant(1.0);
    ActivationStateModel_t state_activation(ps, state_residual.get_nr_dim(), state_activation_weight);
    CostStateModel_t state_cost(ps, state_residual, state_activation);
    cost_manager.addItem("state_regularization", state_cost, 1e-1);

    // Swing foot costs (if any swing feet specified)
    for (size_t i = 0; i < swing_foot_ids.size(); ++i) {
        if (i < swing_foot_targets.size()) {
            ResidualFrameTranslationModel_t foot_residual(ps, swing_foot_ids[i], swing_foot_targets[i]);
            Eigen::Vector<VarScalar_, ResidualFrameTranslationModel_t::DimNR_t::Value> foot_activation_weight =
                Eigen::Vector<VarScalar_, ResidualFrameTranslationModel_t::DimNR_t::Value>::Constant(1.0);
            ActivationFrameTranslationModel_t foot_activation(ps, foot_residual.get_nr_dim(), foot_activation_weight);
            CostFrameTranslationModel_t foot_cost(ps, foot_residual, foot_activation);
            cost_manager.addItem("swing_foot_" + std::to_string(swing_foot_ids[i]), foot_cost, 1e4);

            // Add velocity cost for swing foot
            Motion_t foot_vel_ref = Motion_t::Zero(); // NEED TO MAKE THIS A REAL REFERENCE TO TRACK A FOOT VELOCITY LIKE IN CROCODDYL
            ResidualFrameVelocityModel_t foot_vel_residual(ps, swing_foot_ids[i], foot_vel_ref, pinocchio::LOCAL_WORLD_ALIGNED);
            Eigen::Vector<VarScalar_, ResidualFrameVelocityModel_t::DimNR_t::Value> foot_vel_activation_weight =
                Eigen::Vector<VarScalar_, ResidualFrameVelocityModel_t::DimNR_t::Value>::Constant(1.0);
            ActivationFrameVelocityModel_t foot_vel_activation(ps, foot_vel_residual.get_nr_dim(), foot_vel_activation_weight);
            CostFrameVelocityModel_t foot_vel_cost(ps, foot_vel_residual, foot_vel_activation);
            cost_manager.addItem("swing_foot_vel_" + std::to_string(swing_foot_ids[i]), foot_vel_cost, 1e2);
        }
    }

    return cost_manager;
}

// Function to create contact manager for given support feet
ContactModelManager_t createContactManager(const PhaseSpec_t &ps, const std::vector<pinocchio::FrameIndex> &support_foot_ids)
{
    ContactModelManager_t contact_manager(ps);

    for (const auto &foot_id : support_foot_ids) {
        ContactModel_t contact(ps, foot_id, pinocchio::LOCAL_WORLD_ALIGNED,
                              Eigen::Vector3d::Zero(), Eigen::Vector2d(0., 50.));
        contact_manager.addItem("contact_" + std::to_string(foot_id), contact);
    }

    return contact_manager;
}

int main()
{
    std::string urdf_path = "/home/quant/research/Galileo/resources/go1/urdf/go1.urdf";

    RobotModel_t model = RobotModel_t();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    State_t state = State_t(model);
    ActuationModel_t actuation = ActuationModel_t(state);

    PhaseSpec_t ps = PhaseSpec_t(state);

    std::cout << "Robot initialized with " << model.nq << " positions and " << model.nv << " velocities" << std::endl;

    // Get foot frame IDs
    pinocchio::FrameIndex lf_foot_id = model.getFrameId("LF_FOOT");
    pinocchio::FrameIndex rf_foot_id = model.getFrameId("RF_FOOT");
    pinocchio::FrameIndex lh_foot_id = model.getFrameId("LH_FOOT");
    pinocchio::FrameIndex rh_foot_id = model.getFrameId("RH_FOOT");

    // Setup default state
    Eigen::VectorXd defaultstate = Eigen::VectorXd::Zero(model.nq + model.nv);
    defaultstate.head(model.nq) = model.referenceConfigurations.at("standing");

    // Compute initial foot positions
    pinocchio::Data rdata(model);
    const auto q0 = defaultstate.head(model.nq);
    pinocchio::forwardKinematics(model, rdata, q0);
    pinocchio::centerOfMass(model, rdata, q0);
    pinocchio::updateFramePlacements(model, rdata);

    const Eigen::Vector3d rf_foot_pos0 = rdata.oMf[rf_foot_id].translation();
    const Eigen::Vector3d rh_foot_pos0 = rdata.oMf[rh_foot_id].translation();
    const Eigen::Vector3d lf_foot_pos0 = rdata.oMf[lf_foot_id].translation();
    const Eigen::Vector3d lh_foot_pos0 = rdata.oMf[lh_foot_id].translation();

    Eigen::Vector3d comRef = (rf_foot_pos0 + rh_foot_pos0 + lf_foot_pos0 + lh_foot_pos0) / 4;
    comRef[2] = rdata.com[0][2];

    // Setup control parameters
    JacobiRoots_t jacobi_roots(1.0, 0.0);
    jacobi_roots.compute_roots();
    Eigen::VectorXd nodes = jacobi_roots.get_roots();
    BarycentricInterpolator_t interpolator(nodes);
    ControlParamModel_t control_param(ps, interpolator);

    // Gait parameters
    double timestep = 0.02;
    // TODO: Use these parameters to create multiple phases for a complete walking gait
    // std::size_t stepknots = 20;
    // std::size_t supportknots = 10;
    // double steplength = 0.2;
    // double stepheight = 0.05;

    // Support foot configurations for each swing phase
    std::vector<pinocchio::FrameIndex> rh_support = {lf_foot_id, rf_foot_id, lh_foot_id};
    std::vector<pinocchio::FrameIndex> rf_support = {lf_foot_id, lh_foot_id, rh_foot_id};
    std::vector<pinocchio::FrameIndex> lh_support = {lf_foot_id, rf_foot_id, rh_foot_id};
    std::vector<pinocchio::FrameIndex> lf_support = {rf_foot_id, lh_foot_id, rh_foot_id};

    std::vector<pinocchio::FrameIndex> all_feet_support = {lf_foot_id, rf_foot_id, lh_foot_id, rh_foot_id};

    // Create empty jump model for phases
    ImpulseModelManager_t impulse_manager(ps);
    ConstraintModelManager_t empty_constraint_manager(ps);
    CostModelManager_t empty_cost_manager(ps);
    JumpModel_t jump_model(ps, empty_cost_manager, empty_constraint_manager, impulse_manager);

    // Create phases vector
    std::vector<PhaseModelGeneric_t> phases;

    // Create first double support phase
    auto double_support_cost_manager = createSwingFootCostManager(ps, comRef);
    auto double_support_contact_manager = createContactManager(ps, all_feet_support);

    NodeModel_t double_support_node(ps, double_support_cost_manager, empty_constraint_manager,
                                   double_support_contact_manager, actuation, 0.0, false);

    SegmentModel_t double_support_segment(ps, double_support_node, control_param, timestep);

    PhaseModel_t double_support_phase(ps, jump_model);
    double_support_phase.addSegment(double_support_segment);

    // Add to phases (need to convert to type-erased PhaseModelGeneric_t)
    phases.emplace_back(static_cast<PhaseModelGeneric_t>(double_support_phase));

    std::cout << "Walking gait created successfully!" << std::endl;
    std::cout << "- Robot model has " << model.nq << " DOF" << std::endl;
    std::cout << "- Phase spec created with all required managers" << std::endl;
    std::cout << "- Created " << phases.size() << " phases" << std::endl;

    return 0;
}
