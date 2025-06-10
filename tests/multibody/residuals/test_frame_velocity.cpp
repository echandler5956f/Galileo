#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

#include <galileo/fwd.hpp>
#include <galileo/core/basic-spec.hpp>
#include <galileo/multibody/robot-spec.hpp>

#include <galileo/core/states/state-base.hpp>
#include <galileo/multibody/states/multibody.hpp>

#include <galileo/core/actuations/actuation-base.hpp>
#include <galileo/multibody/actuations/floating-base.hpp>

#include <galileo/predictive/phases/phase-spec.hpp>

#include <galileo/multibody/contacts/implementations/contact-3d.hpp>

#include "galileo/multibody/residuals/frame-velocity.hpp"

#include <galileo/core/constraints/constraint-manager.hpp>
#include <galileo/core/costs/cost-manager.hpp>
#include <galileo/multibody/contacts/contact-manager.hpp>

#include <galileo/predictive/nodes/implementations/node-contact-fwddyn.hpp>

#include <galileo/core/controls/implementations/control-param-jpoly.hpp>

#include <galileo/predictive/segments/implementations/segment-erk-euler.hpp>

#include <galileo/core/data/data-collector-default.hpp>

#include "../../models/pendulum.hpp"
#include "../../helpers/catch_eigen_matchers.hpp"
#include "../../helpers/finite_difference.hpp"

using namespace galileo;
using namespace GalileoMatchers;

// Define a simple RobotSpec for the pendulum
using PendulumBasicSpec = BasicSpecTpl<double, double, 0>;
using PendulumRobotSpec = RobotSpecTpl<PendulumBasicSpec, 0, 1, 0, 1, 0, StateMultibodyTpl, ActuationModelFloatingBaseTpl>;

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

template <typename PS_>
struct DummyPhaseTpl;
template <typename PS_>
struct DummyPhaseModelTpl;
template <typename PS_>
struct DummyPhaseDataTpl;

template <typename PS_>
struct traits<DummyPhaseTpl<PS_>>
{
    using PS = PS_;
    using Model_t = DummyPhaseModelTpl<PS>;
    using Data_t = DummyPhaseDataTpl<PS>;
};

// Define a minimal, but complete, PhaseSpec for the test
template <typename Robot>
using PendulumPhaseSpecTpl = PhaseSpecTpl<Robot,
                                          CostManagerDefaultTpl,
                                          ConstraintManagerDefaultTpl,
                                          NodeTpl,
                                          ControlParamTpl,
                                          SegmentTpl,
                                          DummyPhaseTpl>;
using PendulumPhaseSpec = PendulumPhaseSpecTpl<PendulumRobotSpec>;

TEMPLATE_TEST_CASE("Frame Velocity Residual", "[multibody][residual]", PendulumPhaseSpec)
{
    using PhaseSpec = TestType;
    using RobotSpec = typename PhaseSpec::RS;
    using State = typename RobotSpec::State_t;
    using RobotModel = typename RobotSpec::RobotModel_t;
    using RobotData = typename RobotSpec::RobotData_t;
    using VectorNx = typename RobotSpec::VectorNx_t;
    using VectorNu = typename PhaseSpec::VectorNu_t;

    using ResidualModel = ResidualModelFrameVelocityTpl<PhaseSpec>;
    using ResidualData = typename ResidualModel::Data_t;

    // Setup the model and state
    RobotModel model;
    tests::build_pendulum(model);
    RobotData robot_data(model);

    State state(&model);

    // Get the frame ID for the end-effector
    const std::string frame_name = "link1";
    auto frame_id = model.getFrameId(frame_name);

    // Define a reference velocity
    using Motion = pinocchio::MotionTpl<typename RobotSpec::NumScalar>;
    Motion v_ref(Eigen::Vector3d(0, 0, 1.0), Eigen::Vector3d(0, 0, 0));

    // Create the residual model
    ResidualModel residual_model(&model, frame_id, v_ref, pinocchio::LOCAL);

    using ActuationModel_t = typename RobotSpec::ActuationModel_t;
    using ActuationData_t = typename RobotSpec::ActuationData_t;
    ActuationModel_t actuation;
    ActuationData_t actuation_data = actuation.createData();

    using JointData_t = JointDataTpl<PhaseSpec>;
    JointData_t joint_data(RobotSpec::NV);

    // Create a universal data collector
    using DataCollectorDefault_t = DataCollectorDefaultTpl<PhaseSpec>;
    DataCollectorDefault_t data_collector(&robot_data, &actuation_data, &joint_data);

    GIVEN("A random state and control")
    {
        VectorNx x = state.rand();
        VectorNu u = VectorNu::Random();

        ResidualData residual_data = residual_model.createData(&data_collector);

        pinocchio::forwardKinematics(model, robot_data, x.head(model.nq), x.tail(model.nv));
        pinocchio::updateFramePlacements(model, robot_data);

        WHEN("calc is called")
        {
            residual_model.calc(residual_data, x, u);

            THEN("the residual is the difference to the reference velocity")
            {
                Motion v_actual = pinocchio::getFrameVelocity(model, robot_data, frame_id, pinocchio::LOCAL);
                auto v_diff = v_actual - v_ref;
                REQUIRE_THAT(residual_data.R, Approx(v_diff.toVector(), 1e-9));
            }
        }

        WHEN("calcDiff is called")
        {
            residual_model.calc(residual_data, x, u);

            const auto q = x.derived().template head<RobotSpec::NQ>();
            const typename RobotSpec::VectorNq_t q_vec = q;
            const typename RobotSpec::VectorNv_t v_vec = x.derived().template segment<RobotSpec::NV>(RobotSpec::NQ);
            typename RobotSpec::VectorNv_t a_vec;
            a_vec.setZero();
            pinocchio::computeForwardKinematicsDerivatives(model, robot_data, q_vec, v_vec, a_vec);

            residual_model.calcDiff(residual_data, x, u);

            THEN("the Jacobians match the finite difference approximation")
            {
                typename traits<typename ResidualModel::Meta_t>::Rx_t rx_analytical = residual_data.Rx;
                typename traits<typename ResidualModel::Meta_t>::Ru_t ru_analytical = residual_data.Ru;

                auto func_x = [&](const VectorNx &x_pert, typename traits<typename ResidualModel::Meta_t>::R_t &r_out)
                {
                    pinocchio::forwardKinematics(model, robot_data, x_pert.head(model.nq), x_pert.tail(model.nv));
                    pinocchio::updateFramePlacements(model, robot_data);
                    ResidualData temp_data = residual_model.createData(&data_collector);
                    temp_data.robot = &robot_data;
                    residual_model.calc(temp_data, x_pert, u);
                    r_out = temp_data.R;
                };

                typename traits<typename ResidualModel::Meta_t>::Rx_t rx_numerical(residual_model.nr(), state.get_ndx());
                tests::compute_finite_difference_jacobian_manifold<State, decltype(func_x), VectorNx, typename traits<typename ResidualModel::Meta_t>::R_t, decltype(rx_numerical)>(
                    state, func_x, x, 1e-7, rx_numerical);

                REQUIRE_THAT(rx_analytical, Approx(rx_numerical, 1e-5));
                REQUIRE(ru_analytical.isZero(1e-9));
            }
        }
    }
}