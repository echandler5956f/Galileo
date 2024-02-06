#include "traj_test.h"

int main(int argc, char **argv)
{
    double q0[] = {
        0, 0, 1.0627, 0, 0, 0, 1, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    LeggedBody<Scalar> bot(huron_location, num_ees, end_effector_names);

    Data data = Data(bot.model);

    // Create environment Surfaces
    std::shared_ptr<environment::EnvironmentSurfaces> surfaces = std::make_shared<environment::EnvironmentSurfaces>();
    surfaces->push_back(environment::CreateInfiniteGround());
    std::cout << "EnvironmentSurfaces created" << std::endl;

    std::shared_ptr<contact::ContactSequence> contact_sequence = std::make_shared<contact::ContactSequence>(num_ees);
    std::cout << "ContactSequence created" << std::endl;

    legged::contact::RobotEndEffectors ees = bot.getEndEffectors();
    std::cout << "RobotEndEffectors created" << std::endl;

    ADModel cmodel = bot.model.cast<ADScalar>();

    ADData cdata(cmodel);
    std::cout << "ADModel and ADData created" << std::endl;

    auto nq = bot.model.nq;
    auto nv = bot.model.nv;
    std::shared_ptr<opt::LeggedRobotStates> si = std::make_shared<opt::LeggedRobotStates>(nq, nv, ees);
    std::cout << "LeggedRobotStates created" << std::endl;

    contact::ContactMode initial_mode;
    initial_mode.combination_definition = bot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};
    initial_mode.createModeDynamics(bot.model, ees, si);
    std::cout << "Initial mode created" << std::endl;

    contact_sequence->addPhase(initial_mode, 5, 1.0);
    std::cout << "Initial mode added to sequence" << std::endl;

    casadi::Function F;
    initial_mode.getModeDynamics(F);
    std::cout << "Initial mode dynamics created" << std::endl;

    casadi::SX cx = casadi::SX::sym("x", si->nx);
    casadi::SX cx2 = casadi::SX::sym("x2", si->nx);
    casadi::SX cdx = casadi::SX::sym("dx", si->ndx);
    casadi::SX cu = casadi::SX::sym("u", si->nu);
    casadi::SX cdt = casadi::SX::sym("dt");

    auto ch = si->get_ch(cx);
    auto ch_d = si->get_ch_d(cdx);
    auto cdh = si->get_cdh(cx);
    auto cdh_d = si->get_cdh_d(cdx);
    auto cq = si->get_q(cx);
    auto cq_d = si->get_q_d(cdx);
    auto cqj = si->get_qj(cx);
    auto cv = si->get_v(cx);
    auto cv_d = si->get_v_d(cdx);
    auto cvj = si->get_vj(cx);
    auto cvju = si->get_vju(cu);
    auto cwrenches = si->get_all_wrenches(cu);

    auto ch2 = si->get_ch(cx2);
    auto cdh2 = si->get_cdh(cx2);
    auto cq2 = si->get_q(cx2);
    auto cv2 = si->get_v(cx2);

    ConfigVectorAD q_AD(nq);
    q_AD = Eigen::Map<ConfigVectorAD>(static_cast<std::vector<ADScalar>>(cq).data(), nq, 1);

    ConfigVectorAD q2_AD(nq);
    q2_AD = Eigen::Map<ConfigVectorAD>(static_cast<std::vector<ADScalar>>(cq2).data(), nq, 1);

    opt::TangentVectorAD v_AD(nv);
    v_AD = Eigen::Map<opt::TangentVectorAD>(static_cast<std::vector<opt::ADScalar>>(cv).data(), nv, 1);

    pinocchio::forwardKinematics(cmodel, cdata, q_AD, v_AD); 
    pinocchio::updateFramePlacements(cmodel, cdata);

    casadi::Function Fint("Fint",
                          {cx, cdx, cdt},
                          {vertcat(ch + ch_d,
                                   cdh + cdh_d,
                                   // cq_result, // returns different expression than python version, NO idea why
                                   custom_fint(cx, cdx, cdt),
                                   cv + cv_d)});

    ConfigVectorAD v_result = pinocchio::difference(cmodel, q_AD, q2_AD);
    casadi::SX cv_result(nv, 1);
    pinocchio::casadi::copy(v_result, cv_result);

    casadi::Function Fdif("Fdif",
                          {cx, cx2, cdt},
                          {vertcat(ch2 - ch,
                                   cdh2 - cdh,
                                   cv_result / cdt,
                                   cv2 - cv)});
    casadi::SX cq0(nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    casadi::Function L("L",
                       {cx, cu},
                       {1e-3 * casadi::SX::sumsqr(cvju) +
                        1e-4 * casadi::SX::sumsqr(cwrenches) +
                        1e1 * casadi::SX::sumsqr(cq(casadi::Slice(0, 3)) - cq0(casadi::Slice(0, 3))) +
                        1e1 * casadi::SX::sumsqr(cqj - cq0(casadi::Slice(7, nq)))});

    casadi::Function Phi("Phi",
                         {cx},
                         {1e6 * casadi::SX::sumsqr(cq(casadi::Slice(0, 3)) - cq0(casadi::Slice(0, 3))) +
                          1e6 * casadi::SX::sumsqr(cqj - cq0(casadi::Slice(7, nq)))});

    casadi::Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 100;

    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(Fint, Fdif, F, L, Phi);

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> friction_cone_constraint_builder =
        std::make_shared<FrictionConeConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> velocity_constraint_builder =
        std::make_shared<VelocityConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> contact_constraint_builder =
        std::make_shared<ContactConstraintBuilder<LeggedRobotProblemData>>();

    std::vector<std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>>> builders = {friction_cone_constraint_builder, velocity_constraint_builder};

    std::shared_ptr<LeggedRobotProblemData> legged_problem_data = std::make_shared<LeggedRobotProblemData>(gp_data, surfaces, contact_sequence, si, std::make_shared<ADModel>(cmodel),
                                                                                                           std::make_shared<ADData>(cdata), ees, cx, cu, cdt, 5);

    std::vector<opt::ConstraintData> constraint_datas;
    for (auto builder : builders)
    {
        std::cout << "Building constraint" << std::endl;
        opt::ConstraintData some_data;
        builder->BuildConstraint(*legged_problem_data, 0, some_data);
        std::cout << "Built constraint" << std::endl;
        constraint_datas.push_back(some_data);
    }

    TrajectoryOpt<LeggedRobotProblemData> traj(legged_problem_data, builders, opts);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    traj.init_finite_elements(1, X0);

    casadi::MXVector sol = traj.optimize();

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(100, 0, 1.0);
    Eigen::MatrixXd new_sol = traj.get_solution(new_times);
    Eigen::MatrixXd subMatrix = new_sol.block(si->nh + si->ndh, 0, nq, new_sol.cols());

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