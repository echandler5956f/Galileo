#include "atlas_test.h"

int main(int argc, char **argv)
{
    LeggedBody robot(atlas_location, num_ees, end_effector_names);
    for (int i = 0; i < robot.model.njoints; ++i)
    {
    }
    std::shared_ptr<LeggedRobotStates> si = robot.si;
    ConfigVector q0_vec(robot.model.nq);
    pinocchio::neutral(robot.model, q0_vec);
    q0_vec(2) = 0.95;
    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();
    surfaces->push_back(environment::createInfiniteGround());

    contact::ContactMode initial_mode;
    initial_mode.combination_definition = robot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};

    robot.contact_sequence->addPhase(initial_mode, 20, 1.);

    contact::ContactMode second_mode;
    second_mode.combination_definition = robot.getContactCombination(0b10);
    second_mode.contact_surfaces.resize(second_mode.combination_definition.size());

    // Set the surfaces for the end effectors that are in contact
    int i = 0;
    for (auto combo : second_mode.combination_definition)
    {
        if (combo.second) // If the end effector is in contact
        {
            // Set the surface for this end effector
            second_mode.contact_surfaces[i] = 0;
        }
        else
        {
            second_mode.contact_surfaces[i] = environment::NO_SURFACE;
        }
    }

    // robot.contact_sequence->addPhase(second_mode, 20, 0.1);

    robot.fillModeDynamics();

    casadi::SX cq0(robot.model.nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    Eigen::VectorXd Q_diag(si->ndx);
    Q_diag << 150., 150., 300., 50., 100., 100., /*Centroidal momentum error weights*/
        0., 0., 0., 0., 0., 0.,            /*Rate of Centroidal momentum error weights*/
        500., 500., 500., 0.1, 0.1, 0.1;   /*Floating base position and orientation (exponential coordinates) error weights*/
    for (int i = si->nh + si->ndh + si->nqb; i < si->nh + si->ndh + si->nq; ++i)
        Q_diag(i) = 20.; /*Joint position error weights*/
    for (int i = si->nh + si->ndh + si->nq; i < si->nh + si->ndh + si->nq + si->nvb; ++i)
        Q_diag(i) = 0.; /*Floating base velocity error weights*/
    for (int i = si->nh + si->ndh + si->nq + si->nvb; i < si->nh + si->ndh + si->nq + si->nv; ++i)
        Q_diag(i) = 0.; /*Joint velocity error weights*/
    Eigen::MatrixXd Q_mat = Q_diag.asDiagonal();

    Eigen::VectorXd R_diag(si->nu);
    for (int i = 0; i < si->nF; ++i)
        R_diag(i) = 1e-3; /* Contact wrench error weights*/
    for (int i = si->nF; i < si->nu; ++i)
        R_diag(i) = 10.; /*Joint acceleration error weights*/
    Eigen::MatrixXd R_mat = R_diag.asDiagonal();
    casadi::SX Q = casadi::SX::zeros(si->ndx, si->ndx);
    casadi::SX R = casadi::SX::zeros(si->nu, si->nu);

    pinocchio::casadi::copy(Q_mat, Q);
    pinocchio::casadi::copy(R_mat, R);

    casadi::SX target_pos = vertcat(casadi::SXVector{q0_vec[0] + 0., q0_vec[1] + 0., q0_vec[2] + 0.});
    casadi::SX target_rot = casadi::SX::eye(3);

    pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> oMf = robot.cdata.oMf[robot.model.getFrameId("pelvis", pinocchio::BODY)];
    auto rot = oMf.rotation();
    auto pos = oMf.translation();
    casadi::SX crot = casadi::SX::zeros(3, 3);
    casadi::SX cpos = casadi::SX::zeros(3, 1);
    pinocchio::casadi::copy(rot, crot);
    pinocchio::casadi::copy(pos, cpos);

    casadi::SX rot_c = casadi::SX::mtimes(crot.T(), target_rot);

    casadi::SX target_error_casadi = vertcat(casadi::SXVector{cpos - target_pos, casadi::SX::inv_skew(rot_c - rot_c.T()) / 2});

    casadi::SX X_ref = casadi::SX(X0);
    X_ref(casadi::Slice(si->nh + si->ndh, si->nh + si->ndh + 3)) = target_pos;
    casadi::SX X_error = vertcat(casadi::SXVector{robot.cx(casadi::Slice(0, si->nh + si->ndh)) - X_ref(casadi::Slice(0, si->nh + si->ndh)),
                                                  target_error_casadi,
                                                  robot.cx(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx)) - X_ref(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx))});
    pinocchio::computeTotalMass(robot.model, robot.data);

    casadi::SX U_ref = casadi::SX::zeros(si->nu, 1);
    U_ref(5) = 9.81 * robot.data.mass[0] / 2;
    U_ref(11) = 9.81 * robot.data.mass[0] / 2;
    casadi::SX u_error = robot.cu - U_ref;

    casadi::Function L("L",
                       {robot.cx, robot.cu},
                       {1. * 0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                        1. * 0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});

    casadi::Function Phi("Phi",
                         {robot.cx},
                         {1. * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});

    casadi::Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 250;
    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(robot.fint, robot.fdif, L, Phi);
    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> friction_cone_constraint_builder =
        std::make_shared<FrictionConeConstraintBuilder<LeggedRobotProblemData>>();
    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> velocity_constraint_builder =
        std::make_shared<VelocityConstraintBuilder<LeggedRobotProblemData>>();
    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> contact_constraint_builder =
        std::make_shared<ContactConstraintBuilder<LeggedRobotProblemData>>();
    std::vector<std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>>> builders = {velocity_constraint_builder, friction_cone_constraint_builder};

    std::shared_ptr<LeggedRobotProblemData> legged_problem_data = std::make_shared<LeggedRobotProblemData>(gp_data, surfaces, robot.contact_sequence, si, std::make_shared<ADModel>(robot.cmodel),
                                                                                                           std::make_shared<ADData>(robot.cdata), robot.getEndEffectors(), robot.cx, robot.cu, robot.cdt);

    TrajectoryOpt<LeggedRobotProblemData, contact::ContactMode> traj(legged_problem_data, robot.contact_sequence, builders, opts);

    traj.initFiniteElements(1, X0);
    casadi::MXVector sol = traj.optimize();

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(50, 0, 1.);
    solution_t new_sol = solution_t(new_times);
    traj.getSolution(new_sol);
    auto cons = traj.getConstraintViolations(new_sol);

    for (auto con : cons)
    {
        for (auto c : con)
        {
            std::cout << c.name << std::endl;
            std::cout << c.evaluation_and_bounds << std::endl;
        }
    }

    Eigen::MatrixXd subMatrix = new_sol.state_result.block(si->nh + si->ndh, 0, si->nq, new_sol.state_result.cols());
    // std::cout << subMatrix << std::endl;
    std::ofstream new_times_file("../examples/visualization/sol_times.csv");
    if (new_times_file.is_open())
    {
        new_times_file << new_times.transpose().format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
        new_times_file.close();
    }

    // Save new_sol to a CSV file
    std::ofstream new_sol_states_file("../examples/visualization/sol_states.csv");
    if (new_sol_states_file.is_open())
    {
        new_sol_states_file << subMatrix.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
        new_sol_states_file.close();
    }

    std::ofstream file("../examples/visualization/metadata.csv");
    if (file.is_open())
    {
        file << "urdf location: " << atlas_location << "\n";

        file << "q0: ";
        for (int i = 0; i < q0_vec.size(); ++i)
        {
            file << q0_vec[i];
            if (i != q0_vec.size() - 1) // not the last element
            {
                file << ", ";
            }
        }
        file << "\n";

        file.close();
    }

    return 0;
}