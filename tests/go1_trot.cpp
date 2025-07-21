#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include "galileo/predictive/phases/phase-spec.hpp"

#include "galileo/multibody/actuations/implementations/actuation-floating-base.hpp"
#include "galileo/multibody/states/implementations/state-multibody.hpp"

#include "galileo/core/activations/implementations/activation-quadratic.hpp"
#include "galileo/multibody/residuals/implementations/residual-frame-translation.hpp"

#include "galileo/core/costs/cost-manager.hpp"
#include "galileo/core/costs/fwd.hpp"
#include "galileo/core/costs/implementations/cost-residual.hpp"

#include "galileo/core/constraints/equality/constraint-manager.hpp"
#include "galileo/core/constraints/equality/fwd.hpp"
#include "galileo/core/constraints/equality/implementations/constraint-residual.hpp"

#include "galileo/multibody/contacts/contact-manager.hpp"
#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/contacts/implementations/contact-3d.hpp"

#include "galileo/predictive/nodes/implementations/node-contact-fwddyn.hpp"

#include "galileo/common/math/barycentric-interpolator.hpp"
#include "galileo/common/math/jacobi-roots.hpp"

#include "galileo/core/controls/implementations/control-param-polynomial.hpp"

#include "galileo/predictive/segments/implementations/segment-erk-euler.hpp"

#include "galileo/core/data/data-collector-default.hpp"

#include <iostream>
#include <string>

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
using ConstraintManagerTestTpl = galileo::ConstraintManagerTpl<PhaseSpec, ConstraintCollectionTestTpl>;

template <typename PhaseSpec>
using CostManagerTestTpl = galileo::CostManagerTpl<PhaseSpec, CostCollectionTestTpl>;

template <typename PhaseSpec>
using ContactManagerTestTpl = galileo::ContactManagerTpl<PhaseSpec, ContactCollectionTestTpl>;
template <typename PhaseSpec>
using ContactModelManagerTestTpl = galileo::ContactModelManagerTpl<PhaseSpec, ContactCollectionTestTpl>;

template <typename PhaseSpec>
using NodeTpl = galileo::NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTestTpl>;

static constexpr int NOrder = 1;
template <typename PhaseSpec>
using ControlParamTpl = galileo::ControlParamPolynomialTpl<PhaseSpec, NOrder>;

template <typename PhaseSpec>
using SegmentTpl = galileo::SegmentERKEulerTpl<PhaseSpec>;

template <typename PhaseSpec>
struct DummyPhaseTpl;
template <typename PhaseSpec>
struct DummyPhaseModelTpl;
template <typename PhaseSpec>
struct DummyPhaseDataTpl;

namespace galileo
{
    template <typename PhaseSpec>
    struct traits<DummyPhaseTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;
        using Meta_t = DummyPhaseTpl<PS>;
        using Model_t = DummyPhaseModelTpl<PS>;
        using Data_t = DummyPhaseDataTpl<PS>;
    }; // traits<DummyPhaseTpl<PhaseSpec>>
} // namespace galileo

using PhaseSpec_t = galileo::PhaseSpecTpl<RobotSpec_t, ConstraintManagerTestTpl, CostManagerTestTpl, NodeTpl, ControlParamTpl, SegmentTpl, DummyPhaseTpl>;

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

using NodeModel_t = typename PhaseSpec_t::NodeModel_t;
using NodeData_t = typename PhaseSpec_t::NodeData_t;

using JacobiRoots_t = galileo::JacobiRootsTpl<VarScalar, NOrder, Options>;
using BarycentricInterpolator_t = galileo::BarycentricInterpolatorTpl<VarScalar, NOrder, Options>;

using ControlParamModel_t = typename PhaseSpec_t::ControlParamModel_t;

using SegmentModel_t = typename PhaseSpec_t::SegmentModel_t;
using SegmentData_t = typename PhaseSpec_t::SegmentData_t;

int main(int argc, char *argv[])
{
    std::string urdf_path = "/home/quant/research/Galileo/resources/go1/urdf/go1.urdf";

    RobotModel_t model = RobotModel_t();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);

    State_t state = State_t(model);
    ActuationModel_t actuation = ActuationModel_t(state);

    PhaseSpec_t ps = PhaseSpec_t(state);

    pinocchio::FrameIndex frame_id = model.getFrameId("base");
    ResidualModel_t residual = ResidualModel_t(ps, frame_id, Eigen::Vector3d(0., 0., 0.));
    ActivationModel_t activation = ActivationModel_t(ps, residual.get_nr_dim());

    CostModel_t cost = CostModel_t(ps, residual, activation);
    CostModelManager_t cost_manager = CostModelManager_t(ps);
    cost_manager.addCost("test_cost", cost, 1.0);

    ConstraintModel_t constraint = ConstraintModel_t(ps, residual);
    ConstraintModelManager_t constraint_manager = ConstraintModelManager_t(ps);
    constraint_manager.addConstraint("test_constraint", constraint);

    pinocchio::FrameIndex frame_id_2 = model.getFrameId("LF_FOOT");
    ContactModel_t contact = ContactModel_t(ps, frame_id_2, pinocchio::LOCAL, Eigen::Vector3d(0., 0., 0.), Eigen::Vector2d(0., 0.));
    ContactModelManager_t contact_manager = ContactModelManager_t(ps);
    contact_manager.addContact("test_contact", contact);

    NodeModel_t node = NodeModel_t(ps, cost_manager, constraint_manager, contact_manager, actuation, 0.0, false);
    // NodeData_t node_data = node.createData();

    JacobiRoots_t jacobi_roots = JacobiRoots_t(1.0, 0.0);
    jacobi_roots.compute_roots();
    Eigen::VectorXd nodes = jacobi_roots.get_roots();

    BarycentricInterpolator_t interpolator = BarycentricInterpolator_t(nodes);
    ControlParamModel_t control_param = ControlParamModel_t(ps, interpolator);

    NumScalar period = 0.0;
    SegmentModel_t segment = SegmentModel_t(ps, node, control_param, period);

    SegmentData_t segment_data = segment.createData();

    Eigen::VectorXd x = state.rand();
    Eigen::VectorXd w = Eigen::VectorXd::Zero(ps.get_nw());
    segment.calc(segment_data, x, w);

    return 0;
}
