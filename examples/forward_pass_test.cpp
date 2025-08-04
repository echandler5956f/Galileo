#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "candlewick/multibody/RobotLoader.h"
#include "candlewick/multibody/Visualizer.h"

#include <chrono>
#include <cmath>
#include <algorithm>
#include <pinocchio/algorithm/geometry.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include "galileo/predictive/phases/phase-spec.hpp"

#include "galileo/core/actuations/implementations/actuation-floating-base.hpp"

#include "galileo/core/states/implementations/state-multibody.hpp"

#include "galileo/core/activations/implementations/activation-quadratic.hpp"

#include "galileo/core/residuals/implementations/residual-frame-velocity.hpp"

#include "galileo/core/costs/cost-manager.hpp"
#include "galileo/core/costs/implementations/cost-residual.hpp"

#include "galileo/core/constraints/equality/constraint-manager.hpp"
#include "galileo/core/constraints/equality/implementations/constraint-residual.hpp"

#include "galileo/multibody/contacts/contact-manager.hpp"
#include "galileo/multibody/contacts/implementations/contact-3d.hpp"

#include "galileo/multibody/impulses/impulse-manager.hpp"

#include "galileo/multibody/impulses/implementations/impulse-3d.hpp"

#include "galileo/predictive/nodes/implementations/node-contact-fwddyn.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include "galileo/core/controls/implementations/control-param-polynomial.hpp"

#include "galileo/predictive/segments/implementations/segment-erk-euler.hpp"

#include "galileo/predictive/jumps/implementations/jump-impulse-fwddyn.hpp"

#include "galileo/predictive/phases/implementations/phase-default.hpp"

#include "galileo/core/data/data-collector-default.hpp"

#include "galileo/predictive/phases/phase-generic.hpp"

#include <iostream>
#include <string>

using namespace candlewick::multibody;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

static const RobotSpec go1_robot_spec =
    RobotSpec{
        "urdf/go1.urdf",
        "srdf/go1.srdf",
        fs::path(EXAMPLE_ROBOT_DATA_MODEL_DIR).parent_path(),
        "robots/go1_description",
        true}
        .ensure_absolute_filepaths();

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

template <typename PhaseSpec>
using ResidualFrameVelocityTpl = galileo::ResidualFrameVelocityTpl<PhaseSpec>;

template <typename PhaseSpec>
using ActivationFrameVelocityTpl = galileo::ActivationQuadraticTpl<PhaseSpec, ResidualFrameVelocityTpl>;

template <typename PhaseSpec>
using CostFrameVelocityTpl = galileo::CostResidualTpl<PhaseSpec, ResidualFrameVelocityTpl, ActivationFrameVelocityTpl>;

template <typename PhaseSpec>
using CostFrameVelocityModelTpl = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using CostFrameVelocityDataTpl = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct CostCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using CostModelVariant_t = boost::variant<CostFrameVelocityModelTpl<PS>>;
    using CostDataVariant_t = boost::variant<CostFrameVelocityDataTpl<PS>>;
}; // struct CostCollectionWalkingTpl

template <typename PhaseSpec>
using ConstraintFrameVelocityTpl = galileo::ConstraintResidualTpl<PhaseSpec, ResidualFrameVelocityTpl>;
template <typename PhaseSpec>
using ConstraintModelFrameVelocityTpl = typename galileo::traits<ConstraintFrameVelocityTpl<PhaseSpec>>::Model_t;
template <typename PhaseSpec>
using ConstraintDataFrameVelocityTpl = typename galileo::traits<ConstraintFrameVelocityTpl<PhaseSpec>>::Data_t;

template <typename PhaseSpec>
struct ConstraintCollectionWalkingTpl
{
    using PS = PhaseSpec;
    using ConstraintModelVariant_t = boost::variant<ConstraintModelFrameVelocityTpl<PS>>;
    using ConstraintDataVariant_t = boost::variant<ConstraintDataFrameVelocityTpl<PS>>;
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
using JumpData_t = typename galileo::traits<JumpWalkingTpl<PhaseSpec_t>>::Data_t;
using PhaseModel_t = typename PhaseSpec_t::PhaseModel_t;

using ResidualFrameVelocityModel_t = typename galileo::traits<ResidualFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using ActivationFrameVelocityModel_t = typename galileo::traits<ActivationFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using CostFrameVelocityModel_t = typename galileo::traits<CostFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using ConstraintFrameVelocityModel_t = typename galileo::traits<ConstraintFrameVelocityTpl<PhaseSpec_t>>::Model_t;
using ContactModel_t = typename galileo::traits<ContactWalkingTpl<PhaseSpec_t>>::Model_t;
using ImpulseModel_t = typename galileo::traits<ImpulseWalkingTpl<PhaseSpec_t>>::Model_t;

GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PhaseSpec_t);

int main(int argc, char **argv)
{
    // CLI::App app{"Visualizer example"};
    // argv = app.ensure_utf8(argv);
    // std::array<Uint32, 2> window_dims{1920u, 1080u};
    // double fps;

    // app.add_option("--dims", window_dims, "Window dimensions.")
    //     ->capture_default_str();
    // app.add_option<double, unsigned int>("--fps", fps, "Framerate")
    //     ->default_val(60);

    // CLI11_PARSE(app, argc, argv);

    pinocchio::Model model;
    pinocchio::GeometryModel geom_model;
    loadModels(go1_robot_spec, model, &geom_model, NULL);
    pinocchio::Data rdata(model);

    // Visualizer visualizer{{window_dims[0], window_dims[1]}, model, geom_model};
    // assert(!visualizer.hasExternalData());
    // pinocchio::Data &vis_data = visualizer.data();

    State_t state = State_t(model);

    ActuationModel_t actuation = ActuationModel_t(state);

    PhaseSpec_t ps = PhaseSpec_t(state);

    pinocchio::FrameIndex fr_foot_id = model.getFrameId("FR_foot");
    std::cout << "DEBUG: FR_FOOT found with ID: " << fr_foot_id << std::endl;
    pinocchio::FrameIndex fl_foot_id = model.getFrameId("FL_foot");
    std::cout << "DEBUG: FL_FOOT found with ID: " << fl_foot_id << std::endl;
    pinocchio::FrameIndex rr_foot_id = model.getFrameId("RR_foot");
    std::cout << "DEBUG: RR_FOOT found with ID: " << rr_foot_id << std::endl;
    pinocchio::FrameIndex rl_foot_id = model.getFrameId("RL_foot");
    std::cout << "DEBUG: RL_FOOT found with ID: " << rl_foot_id << std::endl;

    VectorNq_t q0 = model.referenceConfigurations["standing"];
    std::cout << "DEBUG: Standing configuration: " << q0.transpose() << std::endl;
    VectorNv_t v0 = VectorNv_t::Zero(model.nv);
    VectorNx_t x0 = VectorNx_t::Zero(ps.get_nx());
    head(x0, ps.get_nq_dim()) = q0;
    tail(x0, ps.get_nv_dim()) = v0;
    std::cout << "DEBUG: Initial state: " << x0.transpose() << std::endl;

    // Compute initial foot positions with debug output
    pinocchio::forwardKinematics(model, rdata, q0);
    pinocchio::centerOfMass(model, rdata, q0);
    pinocchio::updateFramePlacements(model, rdata);

    const Eigen::Vector3d fr_foot_pos0 = rdata.oMf[fr_foot_id].translation();
    std::cout << "DEBUG: FR_FOOT position: " << fr_foot_pos0.transpose() << std::endl;
    const Eigen::Vector3d fl_foot_pos0 = rdata.oMf[fl_foot_id].translation();
    std::cout << "DEBUG: FL_FOOT position: " << fl_foot_pos0.transpose() << std::endl;
    const Eigen::Vector3d rr_foot_pos0 = rdata.oMf[rr_foot_id].translation();
    std::cout << "DEBUG: RR_FOOT position: " << rr_foot_pos0.transpose() << std::endl;
    const Eigen::Vector3d rl_foot_pos0 = rdata.oMf[rl_foot_id].translation();
    std::cout << "DEBUG: RL_FOOT position: " << rl_foot_pos0.transpose() << std::endl;

    JacobiRoots_t jacobi_roots(1.0, 0.0);
    jacobi_roots.compute_roots();

    Eigen::VectorXd nodes = jacobi_roots.get_roots();
    BarycentricInterpolator_t interpolator(nodes);
    ControlParamModel_t control_param(ps, interpolator);

    // Support foot configurations for each swing phase
    std::vector<pinocchio::FrameIndex> trot_phase_1 = {fr_foot_id, rl_foot_id};
    std::vector<pinocchio::FrameIndex> trot_phase_2 = {fl_foot_id, rr_foot_id};
    std::vector<pinocchio::FrameIndex> all_feet_support = {fr_foot_id, fl_foot_id, rr_foot_id, rl_foot_id};

    ConstraintModelManager_t empty_constraint_manager(ps);
    CostModelManager_t empty_cost_manager(ps);

    std::vector<int> num_knots = {2, 2, 2, 2, 2, 2};
    std::vector<double> phase_durations = {0.5, 0.25, 0.25, 0.25, 0.25, 0.5};
    std::vector<std::vector<pinocchio::FrameIndex>> phases = {all_feet_support, trot_phase_1, trot_phase_2, trot_phase_1, trot_phase_2, all_feet_support};

    std::vector<std::shared_ptr<JumpModel_t>> jump_models;
    std::vector<std::shared_ptr<JumpData_t>> jump_datas;
    std::vector<std::shared_ptr<SegmentModel_t>> segment_models;
    std::vector<std::shared_ptr<SegmentData_t>> segment_datas;

    std::vector<std::shared_ptr<ContactModelManager_t>> contact_managers;
    std::vector<std::shared_ptr<ImpulseModelManager_t>> impulse_managers;

    jump_models.reserve(phases.size());
    jump_datas.reserve(phases.size());
    segment_models.reserve(phases.size() * num_knots.size());
    segment_datas.reserve(phases.size() * num_knots.size());
    contact_managers.reserve(phases.size() * num_knots.size());
    impulse_managers.reserve(phases.size());

    for (int i = 0; i < phases.size(); i++)
    {
        if (i > 0)
        {
            std::cout << "DEBUG: Creating impulse manager for phase " << i << std::endl;
            ImpulseModelManager_t impulse_manager(ps);
            for (const auto &foot_id : phases[i])
            {
                std::cout << "DEBUG: Creating impulse for foot " << foot_id << " phase " << i << std::endl;
                ImpulseModel_t impulse(ps, foot_id, pinocchio::LOCAL_WORLD_ALIGNED);
                impulse_manager.addItem("impulse_" + std::to_string(foot_id), impulse);
                std::cout << "DEBUG: Finished creating impulse for foot " << foot_id << " phase " << i << std::endl;
            }
            impulse_managers.push_back(std::make_shared<ImpulseModelManager_t>(impulse_manager));
        }

        for (int j = 0; j < num_knots[i]; j++)
        {
            std::cout << "DEBUG: Creating contact manager for phase " << i << " knot " << j << std::endl;
            ContactModelManager_t contact_manager(ps);
            for (const auto &foot_id : phases[i])
            {
                std::cout << "DEBUG: Creating contact for foot " << foot_id << " knot " << j << std::endl;
                ContactModel_t contact(ps, foot_id, pinocchio::LOCAL_WORLD_ALIGNED,
                                       Eigen::Vector3d::Zero(), Eigen::Vector2d(0., 50.));
                contact_manager.addItem("contact_" + std::to_string(foot_id), contact);
                std::cout << "DEBUG: Finished creating contact for foot " << foot_id << " knot " << j << std::endl;
            }
            contact_managers.push_back(std::make_shared<ContactModelManager_t>(contact_manager));
        }
    }

    for (int i = 0; i < phases.size(); i++)
    {
        if (i > 0)
        {
            std::cout << "DEBUG: Creating jump model for phase " << i << std::endl;
            JumpModel_t jump_model(ps, empty_cost_manager, empty_constraint_manager, *impulse_managers[i-1]);
            std::cout << "DEBUG: Finished creating jump model for phase " << i << std::endl;
            jump_models.push_back(std::make_shared<JumpModel_t>(jump_model));
            std::cout << "DEBUG: Creating jump data for phase " << i << std::endl;
            jump_datas.push_back(std::make_shared<JumpData_t>(jump_model.createData()));
            std::cout << "DEBUG: Finished creating jump data for phase " << i << std::endl;
        }
        for (int j = 0; j < num_knots[i]; j++)
        {
            std::cout << "DEBUG: Creating node for phase " << i << " knot " << j << std::endl;
            NodeModel_t node(ps, empty_cost_manager, empty_constraint_manager, *contact_managers[i*num_knots[i] + j], actuation, 0.0, false);
            std::cout << "DEBUG: Finished creating node for phase " << i << " knot " << j << std::endl;
            std::cout << "DEBUG: Creating segment model for phase " << i << " knot " << j << std::endl;
            double timestep = phase_durations[i] / num_knots[i];
            SegmentModel_t segment_model(ps, node, control_param, timestep);
            std::cout << "DEBUG: Finished creating segment model for phase " << i << " knot " << j << std::endl;
            segment_models.push_back(std::make_shared<SegmentModel_t>(segment_model));
            std::cout << "DEBUG: Creating segment data for phase " << i << " knot " << j << std::endl;
            segment_datas.push_back(std::make_shared<SegmentData_t>(segment_model.createData()));
            std::cout << "DEBUG: Finished creating segment data for phase " << i << " knot " << j << std::endl;
        }
    }
    std::cout << "DEBUG: Finished creating segment models" << std::endl;

    // Forward simulation pass
    std::vector<VectorNx_t> xs;
    std::vector<VectorNu_t> us;

    VectorNx_t x = VectorNx_t::Zero(state.get_nx());
    head(x, state.get_nq_dim()) = model.referenceConfigurations["standing"];
    tail(x, state.get_nv_dim()) = VectorNv_t::Zero(model.nv);
    std::cout << "DEBUG: Initial state: " << x.transpose() << std::endl;

    std::cout << "DEBUG: Starting forward simulation pass" << std::endl;
    for (int i = 0; i < phases.size(); i++)
    {
        if (i > 0)
        {
            std::cout << "DEBUG: Calculating jump model for phase " << i << std::endl;
            jump_models[i-1]->calc(*jump_datas[i-1], x);
            std::cout << "DEBUG: Finished calculating jump model for phase " << i << std::endl;
            x = jump_datas[i-1]->XNext;
            xs.push_back(x);
        }
        for (int j = 0; j < num_knots[i]; j++)
        {
            Eigen::VectorXd w = Eigen::VectorXd::Zero(ps.get_nw());
            std::cout << "DEBUG: Calculating segment model for phase " << i << " knot " << j << std::endl;
            segment_models[i * num_knots[i] + j]->calc(*segment_datas[i * num_knots[i] + j], x, w);
            std::cout << "DEBUG: Finished calculating segment model for phase " << i << " knot " << j << std::endl;
            x = segment_datas[i * num_knots[i] + j]->XNext;
            xs.push_back(x);
        }
    }

    std::cout << "DEBUG: xs.size(): " << xs.size() << std::endl;
    std::cout << "DEBUG: us.size(): " << us.size() << std::endl;

    // Calculate total trajectory time and create time vector
    double total_time = 0.0;
    std::vector<double> knot_times;
    knot_times.push_back(0.0);

    for (int i = 0; i < phases.size(); i++)
    {
        if (i > 0)
        {
            // Add time for jump (instantaneous)
            knot_times.push_back(total_time);
        }

        double timestep = phase_durations[i] / num_knots[i];
        for (int j = 0; j < num_knots[i]; j++)
        {
            total_time += timestep;
            knot_times.push_back(total_time);
        }
    }

    std::cout << "DEBUG: Total trajectory time: " << total_time << " seconds" << std::endl;
    std::cout << "DEBUG: Number of knot points: " << knot_times.size() << std::endl;

    // double dt = 1. / static_cast<double>(fps);
    // using duration_t = std::chrono::duration<double>;
    // double t = 0.;
    // VectorNx_t q = defaultstate.head(model.nq);
    // VectorNx_t qn = q;
    // VectorNv_t v = VectorNv_t::Zero(model.nv);

    // while (!visualizer.shouldExit())
    // {
    //     const auto now = steady_clock::now();

    //     // Cycle through the trajectory (loop when reaching the end)
    //     double trajectory_t = fmod(t, total_time);

    //     // Find the appropriate knot points for interpolation
    //     int knot_idx = 0;
    //     for (int i = 0; i < knot_times.size() - 1; i++)
    //     {
    //         if (trajectory_t >= knot_times[i] && trajectory_t < knot_times[i + 1])
    //         {
    //             knot_idx = i;
    //             break;
    //         }
    //     }

    //     if (knot_idx >= xs.size() - 1)
    //     {
    //         knot_idx = xs.size() - 2;
    //     }

    //     // Calculate interpolation parameter
    //     double alpha = 0.0;
    //     if (knot_times[knot_idx + 1] > knot_times[knot_idx])
    //     {
    //         alpha = (trajectory_t - knot_times[knot_idx]) / (knot_times[knot_idx + 1] - knot_times[knot_idx]);
    //     }
    //     alpha = std::max(0.0, std::min(1.0, alpha)); // Clamp to [0,1]

    //     // Extract q and v from the state vectors
    //     VectorNx_t q0 = xs[knot_idx].head(model.nq);
    //     VectorNx_t q1 = xs[knot_idx + 1].head(model.nq);
    //     VectorNv_t v0 = xs[knot_idx].tail(model.nv);
    //     VectorNv_t v1 = xs[knot_idx + 1].tail(model.nv);

    //     // Interpolate configuration using pinocchio's manifold interpolation
    //     pinocchio::interpolate(model, q0, q1, alpha, q);

    //     // Linear interpolation for velocities
    //     v = (1.0 - alpha) * v0 + alpha * v1;

    //     // Update kinematics for visualization
    //     pinocchio::forwardKinematics(model, vis_data, q, v);
    //     pinocchio::updateFramePlacements(model, vis_data);

    //     visualizer.display();
    //     std::this_thread::sleep_until(now + duration_t(dt));

    //     t += dt;
    //     qn = q;
    // }

    return 0;
}
