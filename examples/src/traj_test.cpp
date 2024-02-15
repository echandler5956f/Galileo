#include "traj_test.h"

int main(int argc, char **argv)
{
    double q0[] = {
        0, 0, 1.0627, 0, 0, 0, 1, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    LeggedBody robot(huron_location, num_ees, end_effector_names);
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

    robot.contact_sequence->addPhase(second_mode, 20, 0.1);

    robot.fillModeDynamics();

    casadi::SX cq0(robot.model.nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    /*ETH Weights*/
    Eigen::VectorXd Q_diag(si->ndx);
    Q_diag << 15., 15., 30., 5., 10., 10.,                          /*Centroidal momentum error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Rate of Centroidal momentum error weights*/
        500., 500., 500., 200., 200., 200.,                         /*Floating base position and orientation (exponential coordinates) error weights*/
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., /*Joint position error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Floating base velocity error weights*/
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;             /*Joint velocity error weights*/
    Eigen::MatrixXd Q_mat = Q_diag.asDiagonal();

    Eigen::VectorXd R_diag(si->nu);
    R_diag << 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,                   /*First contact wrench error weights*/
        1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,                         /*Second contact wrench error weights*/
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

    // std::vector<opt::ConstraintData> constraint_datas;
    // for (auto builder : builders)
    // {
    //     opt::ConstraintData some_data;
    //     builder->buildConstraint(*legged_problem_data, 0, some_data);
    //     constraint_datas.push_back(some_data);
    // }

    TrajectoryOpt<LeggedRobotProblemData, contact::ContactMode> traj(legged_problem_data, robot.contact_sequence, builders, opts);

    traj.initFiniteElements(1, X0);

    casadi::MXVector sol = traj.optimize();

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(40, 0, 1.1);
    Eigen::MatrixXd new_sol = traj.getSolution(new_times);

    // opt::ConstraintData fri;
    // friction_cone_constraint_builder->buildConstraint(*legged_problem_data, 0, fri);
    // auto force_function = fri.G;
    // opt::ConstraintData vel;
    // velocity_constraint_builder->buildConstraint(*legged_problem_data, 0, vel);
    // auto vel_function = vel.G;

    // auto solu = sol[1];
    // auto dm_u = casadi::MX::evalf(solu(casadi::Slice(0, solu.rows()), casadi::Slice(0, solu.columns())));
    // std::cout << dm_u.get_elements() << std::endl;

    // for (int k = 0; k < new_sol.cols(); ++k)
    // {
    //     casadi::SX tmp_x(si->nx);
    //     casadi::SX tmp_u = casadi::SX::zeros(si->nu, 1);
    //     pinocchio::casadi::copy(new_sol.block(0, k, si->nx, 1), tmp_x);
    //     std::cout << "State: " << tmp_x << std::endl;
    //     casadi::SX tmp_v = vel_function(casadi::SXVector{tmp_x, tmp_u}).at(0);
    //     std::cout << "Velocity: " << tmp_v << std::endl;
    //     casadi::SX tmp_f = force_function(casadi::SXVector{tmp_x, dm_u(casadi::Slice(0, si->nu), k)}).at(0);
    //     std::cout << "Lower bound force: " << fri.lower_bound(casadi::DMVector{0}) << " Force: " << tmp_f << " Upper bound force: " << fri.upper_bound(casadi::DMVector{0}) << std::endl;
    //     auto dm_f = si->get_f<casadi::DM>(dm_u(casadi::Slice(0, si->nu), k), robot.cmodel.getFrameId("l_foot_v_ft_link", pinocchio::BODY));
    //     std::cout << "Actual force: " << dm_f << std::endl;
    // }

    Eigen::MatrixXd subMatrix = new_sol.block(si->nh + si->ndh, 0, si->nq, new_sol.cols());

    std::ofstream new_times_file("../python/new_times.csv");
    if (new_times_file.is_open())
    {
        new_times_file << new_times.transpose().format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
        new_times_file.close();
    }

    // Save new_sol to a CSV file
    std::ofstream new_sol_file("../python/new_sol.csv");
    if (new_sol_file.is_open())
    {
        new_sol_file << subMatrix.format(Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n"));
        new_sol_file.close();
    }

    return 0;
}