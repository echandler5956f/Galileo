#include <pinocchio/fwd.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include <galileo/fwd.hpp>
#include <galileo/core/basic-spec.hpp>
#include <galileo/core/robot-spec.hpp>

#include <galileo/core/states/state-base.hpp>
#include <galileo/core/states/multibody.hpp>

#include <galileo/core/actuations/actuation-base.hpp>
#include <galileo/core/actuations/floating-base.hpp>

#include <galileo/predictive/phases/phase-spec.hpp>

#include <galileo/multibody/contacts/contact-3d.hpp>

#include <galileo/core/constraints/constraint-manager.hpp>
#include <galileo/core/costs/cost-manager.hpp>
#include <galileo/multibody/contacts/contact-manager.hpp>

#include <galileo/predictive/nodes/node-contact-fwddyn.hpp>

#include <galileo/core/controls/control-param-jpoly.hpp>

#include <galileo/predictive/segments/segment-erk-euler.hpp>

// #include <galileo/predictive/trajectory.hpp>

// #include <galileo/predictive/optimal-control-problem.hpp>

// #include <galileo/predictive/solvers/ddp.hpp>

#include <galileo/core/data/data-collector-default.hpp>

#include <iostream>
#include <string>

using namespace galileo;

using VarScalar = double;
using NumScalar = double;
constexpr int Options = Eigen::ColMajor;

using BasicSpec_t = BasicSpecTpl<VarScalar, NumScalar, Options>;

template <typename RS>
using StateTpl = StateMultibodyTpl<RS>;

template <typename RS>
using ActuationTpl = ActuationFloatingBaseTpl<RS>;

constexpr int NQb = 7;
constexpr int NQj = 12;
constexpr int NVb = 6;
constexpr int NVj = 12;
constexpr int NRotors = 0;

// These template classes basically create a code generator for the optimal control problem by means of mixins
using RobotSpec_t = RobotSpecTpl<BasicSpec_t, NQb, NQj, NVb, NVj, NRotors, StateTpl, ActuationTpl>;

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

template <typename PS_>
struct traits<DummyPhaseTpl<PS_>>
{
    using PS = PS_;
    using Model_t = DummyPhaseModelTpl<PS>;
    using Data_t = DummyPhaseDataTpl<PS>;
};

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
    using RobotModel_t = typename PhaseSpec_t::RobotModel_t;
    using RobotData_t = typename PhaseSpec_t::RobotData_t;
    RobotModel_t model = RobotModel_t();
    pinocchio::urdf::buildModel(urdf_path, pinocchio::JointModelFreeFlyerTpl<VarScalar, Options>(), model);
    RobotData_t data = RobotData_t(model);

    StateTpl<RobotSpec_t> state(&model);

    std::cout << "PS::NQb = " << PhaseSpec_t::NQb << std::endl;
    std::cout << "PS::NQj = " << PhaseSpec_t::NQj << std::endl;
    std::cout << "PS::NVb = " << PhaseSpec_t::NVb << std::endl;
    std::cout << "PS::NVj = " << PhaseSpec_t::NVj << std::endl;
    std::cout << "PS::NRotors = " << PhaseSpec_t::NRotors << std::endl;
    std::cout << "PS::NQ = " << PhaseSpec_t::NQ << std::endl;
    std::cout << "PS::NV = " << PhaseSpec_t::NV << std::endl;
    std::cout << "PS::NX = " << PhaseSpec_t::NX << std::endl;
    std::cout << "PS::NDX = " << PhaseSpec_t::NDX << std::endl;
    std::cout << "PS::NUa = " << PhaseSpec_t::NUa << std::endl;
    std::cout << "PS::NU = " << PhaseSpec_t::NU << std::endl;
    std::cout << "PS::NOrder = " << PhaseSpec_t::NOrder << std::endl;
    std::cout << "PS::NW = " << PhaseSpec_t::NW << std::endl;
    std::cout << "PS::NStages = " << PhaseSpec_t::NStages << std::endl;

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

    using JointData_t = JointDataTpl<PhaseSpec_t>;
    JointData_t joint_data(RobotSpec_t::NV);

    // Create a universal data collector
    using DataCollectorDefault_t = DataCollectorDefaultTpl<PhaseSpec_t>;
    DataCollectorDefault_t data_collector(&data, &actuation_data, &joint_data);
    // Test the contact manager
    using ContactModel_t = ContactModel3dTpl<PhaseSpec_t>;
    using ContactData_t = ContactData3dTpl<PhaseSpec_t>;
    ContactModel_t contact(&state, 2, typename PhaseSpec_t::Vector3_t(0, 0, 0), pinocchio::ReferenceFrame::LOCAL, PhaseSpec_t::NU, typename PhaseSpec_t::Vector2_t(0, 50));
    ContactData_t contact_data = contact.createData(&data_collector);
    std::cout << "contact_data.a0_local = " << contact_data.a0_local << std::endl;
    std::cout << "contact_data.fext = " << contact_data.fext << std::endl;

    // Test the contact manager
    using ContactModelManager_t = ContactModelManagerTpl<PhaseSpec_t, ContactCollectionDefaultTpl>;
    using ContactDataManager_t = ContactDataManagerTpl<PhaseSpec_t, ContactCollectionDefaultTpl>;
    ContactModelManager_t contact_model_manager(&state);
    std::string contact_name1 = "contact1";
    std::string contact_name2 = "contact2";
    contact_model_manager.addContact(contact_name1, contact);
    contact_model_manager.addContact(contact_name2, contact);
    contact_model_manager.changeContactStatus(contact_name1, true);
    contact_model_manager.changeContactStatus(contact_name2, false);
    contact_model_manager.removeContact(contact_name1);
    contact_model_manager.changeContactStatus(contact_name2, true);
    contact_model_manager.addContact(contact_name1, contact);

    ContactDataManager_t contact_data_manager = contact_model_manager.createData(&data_collector);
    std::cout << "contact_data_manager.Jc = " << contact_data_manager.Jc << std::endl;
    std::cout << "contact_data_manager.a0 = " << contact_data_manager.a0 << std::endl;
    std::cout << "contact_data_manager.da0_dx = " << contact_data_manager.da0_dx << std::endl;
    std::cout << "contact_data_manager.dv = " << contact_data_manager.dv << std::endl;
    std::cout << "contact_data_manager.ddv_dx = " << contact_data_manager.ddv_dx << std::endl;

    typename PhaseSpec_t::VectorNx_t x0 = state.zero();
    typename PhaseSpec_t::VectorNx_t x = state.rand();
    x.head(3) = x0.head(3);
    std::cout << "x = " << x << std::endl;
    auto q = x.template head<PhaseSpec_t::NQ>();
    auto v = x.template tail<PhaseSpec_t::NV>();
    typename PhaseSpec_t::VectorNu_t u = PhaseSpec_t::VectorNu_t::Zero();

    pinocchio::computeAllTerms(model, *(data_collector.robot), q, v);
    pinocchio::computeCentroidalMomentum(model, *(data_collector.robot));

    actuation.calc(actuation_data, x, u);

    contact_model_manager.calc(contact_data_manager, x);
    std::cout << "contact_data_manager.Jc = " << contact_data_manager.Jc << std::endl;
    std::cout << "contact_data_manager.a0 = " << contact_data_manager.a0 << std::endl;
    std::cout << "contact_data_manager.da0_dx = " << contact_data_manager.da0_dx << std::endl;
    std::cout << "contact_data_manager.dv = " << contact_data_manager.dv << std::endl;
    std::cout << "contact_data_manager.ddv_dx = " << contact_data_manager.ddv_dx << std::endl;

    pinocchio::forwardDynamics(
        model, *(data_collector.robot), actuation_data.tau,
        contact_data_manager.Jc.topRows(contact_model_manager.nc()), contact_data_manager.a0.head(contact_model_manager.nc()),
        VarScalar(0.));
    auto xout = data_collector.robot->ddq;
    contact_model_manager.updateAcceleration(contact_data_manager, data_collector.robot->ddq);
    contact_model_manager.updateForce(contact_data_manager, data_collector.robot->lambda_c);
    data_collector.joint->a = data_collector.robot->ddq;
    data_collector.joint->tau = u;

    Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic> Kinv;
    Kinv.resize(PhaseSpec_t::NV + contact_model_manager.nc(), PhaseSpec_t::NV + contact_model_manager.nc());
    Kinv.setZero();

    pinocchio::computeRNEADerivatives(model, *(data_collector.robot), q, v, xout,
                                    contact_data_manager.fext);
    contact_model_manager.updateRneaDiff(contact_data_manager, *(data_collector.robot));
    pinocchio::getKKTContactDynamicMatrixInverse(
        model, *(data_collector.robot), contact_data_manager.Jc.topRows(contact_model_manager.nc()), Kinv);

    actuation.calcDiff(actuation_data, x, u);

    contact_model_manager.calcDiff(contact_data_manager, x);
    std::cout << "contact_data_manager.Jc = " << contact_data_manager.Jc << std::endl;
    std::cout << "contact_data_manager.a0 = " << contact_data_manager.a0 << std::endl;
    std::cout << "contact_data_manager.da0_dx = " << contact_data_manager.da0_dx << std::endl;
    std::cout << "contact_data_manager.dv = " << contact_data_manager.dv << std::endl;
    std::cout << "contact_data_manager.ddv_dx = " << contact_data_manager.ddv_dx << std::endl;

    // Create a quadruped walking problem
    // 1. Define the phase specs
    // 2. Define the phase collection
    // 3. Define the trajectory
    // 4. Define the OCP
    // 5. Solve the OCP

    return 0;
}
