#include "go1_test.h"

int main(int argc, char **argv)
{
    double q0[] = {
        0., 0., 0.339, 0., 0., 0., 1.,
        0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    LeggedBody robot(go1_location, num_ees, end_effector_names);
    std::shared_ptr<LeggedRobotStates> si = robot.si;

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
    initial_mode.combination_definition = robot.getContactCombination(0b1111);
    initial_mode.contact_surfaces = {0, 0, 0, 0};

    robot.contact_sequence->addPhase(initial_mode, 20, 1.);

    contact::ContactMode second_mode;
    second_mode.combination_definition = robot.getContactCombination(0b1101);
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

    robot.contact_sequence->addPhase(second_mode, 20, 0.25);

    contact::ContactMode final_mode;
    final_mode.combination_definition = robot.getContactCombination(0b1111);
    final_mode.contact_surfaces = {0, 0, 0, 0};

    robot.contact_sequence->addPhase(final_mode, 20, 0.25);

    std::cout << "Filling dynamics" << std::endl;
    robot.fillModeDynamics();

    casadi::SX cq0(robot.model.nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    /*ETH Weights*/
    Eigen::VectorXd Q_diag(si->ndx);
    Q_diag << 15., 15., 30., 5., 10., 10.,                          /*Centroidal momentum error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Rate of Centroidal momentum error weights*/
        500., 500., 500., 0.1, 0.1, 0.1,                         /*Floating base position and orientation (exponential coordinates) error weights*/
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., /*Joint position error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Floating base velocity error weights*/
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;             /*Joint velocity error weights*/
    Eigen::MatrixXd Q_mat = Q_diag.asDiagonal();

    Eigen::VectorXd R_diag(si->nu);
    R_diag << 1e-3, 1e-3, 1e-3,                                     /*First contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Second contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Third contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Fourth contact wrench error weights*/
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.; /*Joint acceleration error weights*/
    Eigen::MatrixXd R_mat = R_diag.asDiagonal();

    casadi::SX Q = casadi::SX::zeros(si->ndx, si->ndx);
    casadi::SX R = casadi::SX::zeros(si->nu, si->nu);

    pinocchio::casadi::copy(Q_mat, Q);
    pinocchio::casadi::copy(R_mat, R);

    casadi::SX target_pos = vertcat(casadi::SXVector{q0[0] + 0., q0[1] + 0., q0[2] + 0.});
    casadi::SX target_rot = casadi::SX::eye(3);

    pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> oMf = robot.cdata.oMf[robot.model.getFrameId("base", pinocchio::BODY)];
    auto rot = oMf.rotation();
    auto pos = oMf.translation();
    casadi::SX crot = casadi::SX::zeros(3, 3);
    casadi::SX cpos = casadi::SX::zeros(3, 1);
    pinocchio::casadi::copy(rot, crot);
    pinocchio::casadi::copy(pos, cpos);

    casadi::SX rot_c = casadi::SX::mtimes(crot.T(), target_rot);

    casadi::SX target_error_casadi = vertcat(casadi::SXVector{cpos - target_pos, casadi::SX::inv_skew(rot_c - rot_c.T()) / 2});

    std::cout << "Target error: " << target_error_casadi << std::endl;

    casadi::SX X_ref = casadi::SX(X0);
    X_ref(casadi::Slice(si->nh + si->ndh, si->nh + si->ndh + 3)) = target_pos;
    casadi::SX X_error = vertcat(casadi::SXVector{robot.cx(casadi::Slice(0, si->nh + si->ndh)) - X_ref(casadi::Slice(0, si->nh + si->ndh)),
                                                  target_error_casadi,
                                                  robot.cx(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx)) - X_ref(casadi::Slice(si->nh + si->ndh + si->nqb, si->nx))});

    pinocchio::computeTotalMass(robot.model, robot.data);

    casadi::SX U_ref = casadi::SX::zeros(si->nu, 1);
    // U_ref(2) = 9.81 * robot.data.mass[0] / robot.num_end_effectors_;
    // U_ref(5) = 9.81 * robot.data.mass[0] / robot.num_end_effectors_;
    // U_ref(8) = 9.81 * robot.data.mass[0] / robot.num_end_effectors_;
    // U_ref(11) = 9.81 * robot.data.mass[0] / robot.num_end_effectors_;
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

    std::cout << "Total duration: " << robot.contact_sequence->getDT() << std::endl;
    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(50, 0., robot.contact_sequence->getDT());
    solution_t new_sol = solution_t(new_times);
    traj.getSolution(new_sol);

     auto cons = traj.getConstraintViolations(new_sol);

    for (auto con : cons)
    {
        for (auto c : con)
        {
            std::cout << c.name << ": " << std::endl;
            auto evaled = c.evaluation_and_bounds.transpose();
            std::cout << evaled << std::endl;
        }
    }

    Eigen::MatrixXd subMatrix = new_sol.state_result.block(si->nh + si->ndh, 0, si->nq, new_sol.state_result.cols());

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
        file << "urdf location: " << go1_location << "\n";

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