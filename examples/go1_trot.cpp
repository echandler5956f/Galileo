#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "galileo/predictive/phases/phase-spec.hpp"

#include "galileo/domains/multibody/core/actuations/impl/actuation-floating-base.hpp"

#include "galileo/domains/multibody/core/states/impl/state-multibody.hpp"

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
#include "galileo/domains/multibody/spatial/contacts/contact-manager.hpp"

#include "galileo/domains/multibody/spatial/impulses/impl/impulse-3d.hpp"
#include "galileo/domains/multibody/spatial/impulses/impulse-manager.hpp"

#include "galileo/domains/multibody/predictive/nodes/impl/node-contact-fwddyn.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include "galileo/core/controls/impl/control-param-polynomial.hpp"

#include "galileo/predictive/segments/impl/segment-erk-euler.hpp"

#include "galileo/domains/multibody/predictive/jumps/impl/jump-impulse-fwddyn.hpp"

#include "galileo/predictive/phases/impl/phase-default.hpp"

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
    std::cout << "DEBUG: Entering createSwingFootCostManager" << std::endl;
    std::cout << "DEBUG: COM reference: " << com_ref.transpose() << std::endl;
    std::cout << "DEBUG: Number of swing feet: " << swing_foot_ids.size() << std::endl;
    std::cout << "DEBUG: Number of swing targets: " << swing_foot_targets.size() << std::endl;

    std::cout << "DEBUG: Creating cost manager..." << std::endl;
    CostModelManager_t cost_manager(ps);
    std::cout << "DEBUG: Cost manager created" << std::endl;

    // CoM position cost (always active)
    std::cout << "DEBUG: Creating CoM residual..." << std::endl;
    ResidualCoMPositionModel_t com_residual(ps, com_ref);
    std::cout << "DEBUG: CoM residual created" << std::endl;

    std::cout << "DEBUG: Creating CoM activation weight..." << std::endl;
    Eigen::Vector<VarScalar_, ResidualCoMPositionModel_t::DimNR_t::Value> com_activation_weight =
        Eigen::Vector<VarScalar_, ResidualCoMPositionModel_t::DimNR_t::Value>::Constant(1.0);
    std::cout << "DEBUG: CoM activation weight created with size: " << com_activation_weight.size() << std::endl;

    std::cout << "DEBUG: Creating CoM activation..." << std::endl;
    ActivationCoMModel_t com_activation(ps, com_residual.get_nr_dim(), com_activation_weight);
    std::cout << "DEBUG: CoM activation created" << std::endl;

    std::cout << "DEBUG: Creating CoM cost..." << std::endl;
    CostCoMPositionModel_t com_cost(ps, com_residual, com_activation);
    std::cout << "DEBUG: CoM cost created" << std::endl;

    std::cout << "DEBUG: Adding CoM cost to manager..." << std::endl;
    cost_manager.addItem("com_position", com_cost, 1e3);
    std::cout << "DEBUG: CoM cost added" << std::endl;

    // Control regularization cost
    std::cout << "DEBUG: Creating control regularization..." << std::endl;
    Eigen::Vector<VarScalar_, PhaseSpec_t::DimNU_t::Value> control_ref =
        Eigen::Vector<VarScalar_, PhaseSpec_t::DimNU_t::Value>::Zero();
    std::cout << "DEBUG: Control reference size: " << control_ref.size() << std::endl;

    ResidualControlModel_t control_residual(ps, control_ref);
    std::cout << "DEBUG: Control residual created" << std::endl;

    Eigen::Vector<VarScalar_, ResidualControlModel_t::DimNR_t::Value> control_activation_weight =
        Eigen::Vector<VarScalar_, ResidualControlModel_t::DimNR_t::Value>::Constant(1.0);
    ActivationControlModel_t control_activation(ps, control_residual.get_nr_dim(), control_activation_weight);
    CostControlModel_t control_cost(ps, control_residual, control_activation);
    cost_manager.addItem("control_regularization", control_cost, 1e-1);
    std::cout << "DEBUG: Control cost added" << std::endl;

    // State regularization cost
    std::cout << "DEBUG: Creating state regularization..." << std::endl;
    ResidualStateModel_t state_residual(ps, ps.get_state().zero());
    std::cout << "DEBUG: State residual created" << std::endl;

    Eigen::Vector<VarScalar_, ResidualStateModel_t::DimNR_t::Value> state_activation_weight =
        Eigen::Vector<VarScalar_, ResidualStateModel_t::DimNR_t::Value>::Constant(1.0);
    ActivationStateModel_t state_activation(ps, state_residual.get_nr_dim(), state_activation_weight);
    CostStateModel_t state_cost(ps, state_residual, state_activation);
    cost_manager.addItem("state_regularization", state_cost, 1e-1);
    std::cout << "DEBUG: State cost added" << std::endl;

    // Swing foot costs (if any swing feet specified)
    std::cout << "DEBUG: Processing swing foot costs..." << std::endl;
    for (size_t i = 0; i < swing_foot_ids.size(); ++i) {
        std::cout << "DEBUG: Processing swing foot " << i << " with ID: " << swing_foot_ids[i] << std::endl;
        if (i < swing_foot_targets.size()) {
            std::cout << "DEBUG: Target position: " << swing_foot_targets[i].transpose() << std::endl;

            ResidualFrameTranslationModel_t foot_residual(ps, swing_foot_ids[i], swing_foot_targets[i]);
            std::cout << "DEBUG: Foot residual created" << std::endl;

            Eigen::Vector<VarScalar_, ResidualFrameTranslationModel_t::DimNR_t::Value> foot_activation_weight =
                Eigen::Vector<VarScalar_, ResidualFrameTranslationModel_t::DimNR_t::Value>::Constant(1.0);
            ActivationFrameTranslationModel_t foot_activation(ps, foot_residual.get_nr_dim(), foot_activation_weight);
            CostFrameTranslationModel_t foot_cost(ps, foot_residual, foot_activation);
            cost_manager.addItem("swing_foot_" + std::to_string(swing_foot_ids[i]), foot_cost, 1e4);
            std::cout << "DEBUG: Foot position cost added" << std::endl;

            // Add velocity cost for swing foot
            Motion_t foot_vel_ref = Motion_t::Zero(); // NEED TO MAKE THIS A REAL REFERENCE TO TRACK A FOOT VELOCITY LIKE IN CROCODDYL
            std::cout << "DEBUG: Creating foot velocity cost..." << std::endl;
            ResidualFrameVelocityModel_t foot_vel_residual(ps, swing_foot_ids[i], foot_vel_ref, pinocchio::LOCAL_WORLD_ALIGNED);
            std::cout << "DEBUG: Foot velocity residual created" << std::endl;

            Eigen::Vector<VarScalar_, ResidualFrameVelocityModel_t::DimNR_t::Value> foot_vel_activation_weight =
                Eigen::Vector<VarScalar_, ResidualFrameVelocityModel_t::DimNR_t::Value>::Constant(1.0);
            ActivationFrameVelocityModel_t foot_vel_activation(ps, foot_vel_residual.get_nr_dim(), foot_vel_activation_weight);
            CostFrameVelocityModel_t foot_vel_cost(ps, foot_vel_residual, foot_vel_activation);
            cost_manager.addItem("swing_foot_vel_" + std::to_string(swing_foot_ids[i]), foot_vel_cost, 1e2);
            std::cout << "DEBUG: Foot velocity cost added" << std::endl;
        }
    }

    std::cout << "DEBUG: Exiting createSwingFootCostManager" << std::endl;
    return cost_manager;
}

// Function to create contact manager for given support feet
ContactModelManager_t createContactManager(const PhaseSpec_t &ps, const std::vector<pinocchio::FrameIndex> &support_foot_ids)
{
    std::cout << "DEBUG: Entering createContactManager" << std::endl;
    std::cout << "DEBUG: Number of support feet: " << support_foot_ids.size() << std::endl;

    ContactModelManager_t contact_manager(ps);
    std::cout << "DEBUG: Contact manager created" << std::endl;

    for (const auto &foot_id : support_foot_ids) {
        std::cout << "DEBUG: Creating contact for foot ID: " << foot_id << std::endl;
        ContactModel_t contact(ps, foot_id, pinocchio::LOCAL_WORLD_ALIGNED,
                              Eigen::Vector3d::Zero(), Eigen::Vector2d(0., 50.));
        std::cout << "DEBUG: Contact model created" << std::endl;

        contact_manager.addItem("contact_" + std::to_string(foot_id), contact);
        std::cout << "DEBUG: Contact added to manager" << std::endl;
    }

    std::cout << "DEBUG: Exiting createContactManager" << std::endl;
    return contact_manager;
}

int main()
{
    std::cout << "=== DEBUG: Starting go1_trot ===" << std::endl;

    std::string urdf_path = "/home/quant/research/Galileo/resources/go1/urdf/go1.urdf";
    std::cout << "DEBUG: URDF path set to: " << urdf_path << std::endl;

    std::cout << "DEBUG: Creating robot model..." << std::endl;
    RobotModel_t model = RobotModel_t();
    std::cout << "DEBUG: Robot model created successfully" << std::endl;

    std::cout << "DEBUG: Building model from URDF..." << std::endl;
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar_, Options_>(), model);
    std::cout << "DEBUG: Model built successfully" << std::endl;

    std::cout << "DEBUG: Creating state..." << std::endl;
    State_t state = State_t(model);
    std::cout << "DEBUG: State created successfully" << std::endl;

    std::cout << "DEBUG: Creating actuation..." << std::endl;
    ActuationModel_t actuation = ActuationModel_t(state);
    std::cout << "DEBUG: Actuation created successfully" << std::endl;

    std::cout << "DEBUG: Creating phase spec..." << std::endl;
    PhaseSpec_t ps = PhaseSpec_t(state);
    std::cout << "DEBUG: Phase spec created successfully" << std::endl;

    std::cout << "Robot initialized with " << model.nq << " positions and " << model.nv << " velocities" << std::endl;

    // Get foot frame IDs with debug output
    std::cout << "DEBUG: Getting foot frame IDs..." << std::endl;
    std::cout << "DEBUG: Available frames in model:" << std::endl;
    for (size_t i = 0; i < model.frames.size(); ++i) {
        std::cout << "  Frame " << i << ": " << model.frames[i].name << std::endl;
    }

    std::cout << "DEBUG: Looking for LF_FOOT frame..." << std::endl;
    pinocchio::FrameIndex lf_foot_id;
    try {
        lf_foot_id = model.getFrameId("LF_FOOT");
        std::cout << "DEBUG: LF_FOOT found with ID: " << lf_foot_id << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: LF_FOOT not found: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "DEBUG: Looking for RF_FOOT frame..." << std::endl;
    pinocchio::FrameIndex rf_foot_id;
    try {
        rf_foot_id = model.getFrameId("RF_FOOT");
        std::cout << "DEBUG: RF_FOOT found with ID: " << rf_foot_id << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: RF_FOOT not found: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "DEBUG: Looking for LH_FOOT frame..." << std::endl;
    pinocchio::FrameIndex lh_foot_id;
    try {
        lh_foot_id = model.getFrameId("LH_FOOT");
        std::cout << "DEBUG: LH_FOOT found with ID: " << lh_foot_id << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: LH_FOOT not found: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "DEBUG: Looking for RH_FOOT frame..." << std::endl;
    pinocchio::FrameIndex rh_foot_id;
    try {
        rh_foot_id = model.getFrameId("RH_FOOT");
        std::cout << "DEBUG: RH_FOOT found with ID: " << rh_foot_id << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: RH_FOOT not found: " << e.what() << std::endl;
        return -1;
    }

    // Setup default state with debug output
    std::cout << "DEBUG: Setting up default state..." << std::endl;
    std::cout << "DEBUG: Available reference configurations:" << std::endl;
    for (const auto& config : model.referenceConfigurations) {
        std::cout << "  Config: " << config.first << " (size: " << config.second.size() << ")" << std::endl;
    }

    Eigen::VectorXd defaultstate = Eigen::VectorXd::Zero(model.nq + model.nv);
    std::cout << "DEBUG: Default state vector created with size: " << defaultstate.size() << std::endl;

    std::cout << "DEBUG: Looking for 'standing' configuration..." << std::endl;
    try {
        defaultstate.head(model.nq) = model.referenceConfigurations.at("standing");
        std::cout << "DEBUG: Standing configuration loaded successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR: Standing configuration not found: " << e.what() << std::endl;
        std::cout << "DEBUG: Using zero configuration instead" << std::endl;
        defaultstate.head(model.nq) = Eigen::VectorXd::Zero(model.nq);
    }

    // Compute initial foot positions with debug output
    std::cout << "DEBUG: Computing initial foot positions..." << std::endl;
    pinocchio::Data rdata(model);
    const auto q0 = defaultstate.head(model.nq);
    std::cout << "DEBUG: Configuration q0 size: " << q0.size() << std::endl;

    std::cout << "DEBUG: Running forward kinematics..." << std::endl;
    pinocchio::forwardKinematics(model, rdata, q0);
    std::cout << "DEBUG: Forward kinematics completed" << std::endl;

    std::cout << "DEBUG: Computing center of mass..." << std::endl;
    pinocchio::centerOfMass(model, rdata, q0);
    std::cout << "DEBUG: Center of mass computed" << std::endl;

    std::cout << "DEBUG: Updating frame placements..." << std::endl;
    pinocchio::updateFramePlacements(model, rdata);
    std::cout << "DEBUG: Frame placements updated" << std::endl;

    std::cout << "DEBUG: Extracting foot positions..." << std::endl;
    const Eigen::Vector3d rf_foot_pos0 = rdata.oMf[rf_foot_id].translation();
    std::cout << "DEBUG: RF foot position: " << rf_foot_pos0.transpose() << std::endl;

    const Eigen::Vector3d rh_foot_pos0 = rdata.oMf[rh_foot_id].translation();
    std::cout << "DEBUG: RH foot position: " << rh_foot_pos0.transpose() << std::endl;

    const Eigen::Vector3d lf_foot_pos0 = rdata.oMf[lf_foot_id].translation();
    std::cout << "DEBUG: LF foot position: " << lf_foot_pos0.transpose() << std::endl;

    const Eigen::Vector3d lh_foot_pos0 = rdata.oMf[lh_foot_id].translation();
    std::cout << "DEBUG: LH foot position: " << lh_foot_pos0.transpose() << std::endl;

    std::cout << "DEBUG: Computing COM reference..." << std::endl;
    Eigen::Vector3d comRef = (rf_foot_pos0 + rh_foot_pos0 + lf_foot_pos0 + lh_foot_pos0) / 4;
    comRef[2] = rdata.com[0][2];
    std::cout << "DEBUG: COM reference: " << comRef.transpose() << std::endl;

    // Setup control parameters with debug output
    std::cout << "DEBUG: Setting up control parameters..." << std::endl;
    JacobiRoots_t jacobi_roots(1.0, 0.0);
    std::cout << "DEBUG: Jacobi roots created" << std::endl;

    jacobi_roots.compute_roots();
    std::cout << "DEBUG: Jacobi roots computed" << std::endl;

    Eigen::VectorXd nodes = jacobi_roots.get_roots();
    std::cout << "DEBUG: Nodes extracted, size: " << nodes.size() << std::endl;

    BarycentricInterpolator_t interpolator(nodes);
    std::cout << "DEBUG: Barycentric interpolator created" << std::endl;

    ControlParamModel_t control_param(ps, interpolator);
    std::cout << "DEBUG: Control parameter model created" << std::endl;

    // Gait parameters
    double timestep = 0.02;
    std::cout << "DEBUG: Timestep set to: " << timestep << std::endl;

    // Support foot configurations for each swing phase
    std::cout << "DEBUG: Setting up support configurations..." << std::endl;
    std::vector<pinocchio::FrameIndex> rh_support = {lf_foot_id, rf_foot_id, lh_foot_id};
    std::vector<pinocchio::FrameIndex> rf_support = {lf_foot_id, lh_foot_id, rh_foot_id};
    std::vector<pinocchio::FrameIndex> lh_support = {lf_foot_id, rf_foot_id, rh_foot_id};
    std::vector<pinocchio::FrameIndex> lf_support = {rf_foot_id, lh_foot_id, rh_foot_id};
    std::vector<pinocchio::FrameIndex> all_feet_support = {lf_foot_id, rf_foot_id, lh_foot_id, rh_foot_id};
    std::cout << "DEBUG: Support configurations created" << std::endl;

    // Create empty jump model for phases
    std::cout << "DEBUG: Creating managers..." << std::endl;
    ImpulseModelManager_t impulse_manager(ps);
    std::cout << "DEBUG: Impulse manager created" << std::endl;

    ConstraintModelManager_t empty_constraint_manager(ps);
    std::cout << "DEBUG: Constraint manager created" << std::endl;

    CostModelManager_t empty_cost_manager(ps);
    std::cout << "DEBUG: Cost manager created" << std::endl;

    std::cout << "DEBUG: Creating jump model..." << std::endl;
    JumpModel_t jump_model(ps, empty_cost_manager, empty_constraint_manager, impulse_manager);
    std::cout << "DEBUG: Jump model created" << std::endl;

    // Create phases vector
    std::cout << "DEBUG: Creating phases vector..." << std::endl;
    std::vector<PhaseModelGeneric_t> phase_models;
    std::cout << "DEBUG: Phases vector created" << std::endl;

    std::cout << "DEBUG: Creating phase data vector..." << std::endl;
    std::vector<PhaseDataGeneric_t> phase_data;
    std::cout << "DEBUG: Phase data vector created" << std::endl;

    // Create first double support phase
    std::cout << "DEBUG: Creating double support cost manager..." << std::endl;
    auto double_support_cost_manager = createSwingFootCostManager(ps, comRef);
    std::cout << "DEBUG: Double support cost manager created" << std::endl;

    std::cout << "DEBUG: Creating double support contact manager..." << std::endl;
    auto double_support_contact_manager = createContactManager(ps, all_feet_support);
    std::cout << "DEBUG: Double support contact manager created" << std::endl;

    std::cout << "DEBUG: Creating double support node..." << std::endl;
    NodeModel_t double_support_node(ps, double_support_cost_manager, empty_constraint_manager,
                                   double_support_contact_manager, actuation, 0.0, false);
    std::cout << "DEBUG: Double support node created" << std::endl;

    std::cout << "DEBUG: Creating double support segment..." << std::endl;
    SegmentModel_t double_support_segment(ps, double_support_node, control_param, timestep);
    std::cout << "DEBUG: Double support segment created" << std::endl;

    std::cout << "DEBUG: Creating double support phase..." << std::endl;
    PhaseModel_t double_support_phase(ps, jump_model);
    std::cout << "DEBUG: Double support phase created" << std::endl;

    std::cout << "DEBUG: Adding segment to phase..." << std::endl;
    double_support_phase.addSegment(double_support_segment);
    std::cout << "DEBUG: Segment added to phase" << std::endl;

    // Add to phases (need to convert to type-erased PhaseModelGeneric_t)
    std::cout << "DEBUG: Converting and adding phase to vector..." << std::endl;
    phase_models.emplace_back(static_cast<PhaseModelGeneric_t>(double_support_phase));
    std::cout << "DEBUG: Phase added to vector" << std::endl;

    std::cout << "Walking gait created successfully!" << std::endl;
    std::cout << "- Robot model has " << model.nq << " DOF" << std::endl;
    std::cout << "- Phase spec created with all required managers" << std::endl;
    std::cout << "- Created " << phase_models.size() << " phases" << std::endl;

    std::cout << "DEBUG: Creating phase data..." << std::endl;
    for (size_t i = 0; i < phase_models.size(); i++) {
        phase_data.emplace_back(static_cast<PhaseDataGeneric_t>(phase_models[i].createData()));
        std::cout << "DEBUG: Phase data created for phase model at index " << i << std::endl;
    }
    std::cout << "DEBUG: Phase data created" << std::endl;

    std::cout << "DEBUG: Running phase calc..." << std::endl;
    for (size_t i = 0; i < phase_models.size(); i++) {
        Eigen::MatrixXd xs = Eigen::MatrixXd::Zero(ps.get_nx(), 1); // only added one segment
        Eigen::MatrixXd ws = Eigen::MatrixXd::Zero(ps.get_nw(), 1);
        phase_models[i].calc(phase_data[i], xs, ws);
        std::cout << "DEBUG: Phase calc completed for phase model at index " << i << std::endl;
    }
    std::cout << "DEBUG: Finished running phase calc" << std::endl;

    std::cout << "DEBUG: Running phase calc diff..." << std::endl;
    for (size_t i = 0; i < phase_models.size(); i++) {
        Eigen::MatrixXd xs = Eigen::MatrixXd::Zero(ps.get_nx(), 1); // only added one segment
        Eigen::MatrixXd ws = Eigen::MatrixXd::Zero(ps.get_nw(), 1);
        phase_models[i].calcDiff(phase_data[i], xs, ws);
        std::cout << "DEBUG: Phase calc diff completed for phase model at index " << i << std::endl;
    }
    std::cout << "DEBUG: Finished running phase calc diff" << std::endl;

    std::cout << "=== DEBUG: Program completed successfully ===" << std::endl;
    return 0;
}
