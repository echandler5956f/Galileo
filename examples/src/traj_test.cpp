#include "traj_test.h"

int main()
{
    double q0[] = {
        0, 0, 1.0627, 0, 0, 0, 1, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    legged::LeggedBody<Scalar> bot(huron_location, num_ees, end_effector_names);

    Model model = bot.model;
    Data data = Data(model);

    std::cout << "Set the robot" << std::endl;

    // Create environment Surfaces
    std::shared_ptr<legged::environment::EnvironmentSurfaces> surfaces = std::make_shared<legged::environment::EnvironmentSurfaces>();
    surfaces->push_back(legged::environment::CreateInfiniteGround());

    std::cout << "Set the ground" << std::endl;

    std::shared_ptr<galileo::legged::contact::ContactSequence> contact_sequence = make_shared<legged::contact::ContactSequence>(num_ees);

    legged::contact::ContactMode initial_mode;
    initial_mode.combination_definition = bot.getContactCombination(0b11);
    initial_mode.contact_surfaces = {0, 0};

    contact_sequence->addPhase(initial_mode, 100, 0.2);

    pinocchio::computeTotalMass(model, data);
    pinocchio::framesForwardKinematics(model, data, q0_vec);

    ADModel cmodel = model.cast<ADScalar>();
    ADData cdata(cmodel);

    auto mass = data.mass[0];
    auto g = SX::zeros(3, 1);
    g(2) = 9.81;
    auto nq = model.nq;
    auto nv = model.nv;
    shared_ptr<opt::LeggedRobotStates> si = make_shared<opt::LeggedRobotStates>(std::vector<int>{nq, nv});

    SX cx = SX::sym("x", si->nx);
    SX cx2 = SX::sym("x2", si->nx);
    SX cdx = SX::sym("dx", si->ndx);
    SX cu = SX::sym("u", si->nu);
    SX cdt = SX::sym("dt");

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
    auto cf = si->get_f(cu);
    auto ctau = si->get_tau(cu);
    auto cvju = si->get_vju(cu);

    ConfigVectorAD q_AD(model.nq);
    TangentVectorAD v_AD(model.nv);
    TangentVectorAD q_d_AD(model.nv);

    q_AD = Eigen::Map<ConfigVectorAD>(static_cast<vector<ADScalar>>(cq).data(), model.nq, 1);
    v_AD = Eigen::Map<TangentVectorAD>(static_cast<vector<ADScalar>>(cv).data(), model.nv, 1);
    q_d_AD = Eigen::Map<TangentVectorAD>(static_cast<vector<ADScalar>>(cq_d).data(), model.nv, 1);

    pinocchio::centerOfMass(cmodel, cdata, q_AD, false);
    pinocchio::computeCentroidalMap(cmodel, cdata, q_AD);
    pinocchio::forwardKinematics(cmodel, cdata, q_AD, v_AD);
    pinocchio::updateFramePlacements(cmodel, cdata);

    ConfigVectorAD q_result = pinocchio::integrate(cmodel, q_AD, q_d_AD * cdt);
    SX cq_result(model.nq, 1);
    pinocchio::casadi::copy(q_result, cq_result);

    Function Fint("Fint",
                  {cx, cdx, cdt},
                  {vertcat(ch + ch_d,
                           cdh + cdh_d,
                           // cq_result, // returns different expression than python version, NO idea why
                           custom_fint(cx, cdx, cdt),
                           cv + cv_d)});

    auto ch2 = si->get_ch(cx2);
    auto cdh2 = si->get_cdh(cx2);
    auto cq2 = si->get_q(cx2);
    auto cv2 = si->get_v(cx2);

    ConfigVectorAD q2_AD(model.nq);
    q2_AD = Eigen::Map<ConfigVectorAD>(static_cast<vector<ADScalar>>(cq2).data(), model.nq, 1);

    ConfigVectorAD v_result = pinocchio::difference(cmodel, q_AD, q2_AD);
    SX cv_result(model.nv, 1);
    pinocchio::casadi::copy(v_result, cv_result);

    Function Fdif("Fdif",
                  {cx, cx2, cdt},
                  {vertcat(ch2 - ch,
                           cdh2 - cdh,
                           cv_result / cdt,
                           cv2 - cv)});

    auto Ag = cdata.Ag;
    SX cAg(Ag.rows(), Ag.cols());
    pinocchio::casadi::copy(Ag, cAg);

    Function F("F",
               {cx, cu},
               {vertcat(cdh,
                        (cf - mass * g) / mass,
                        ctau / mass,
                        cv,
                        SX::mtimes(SX::inv(cAg(Slice(0, 6), Slice(0, 6))), (mass * ch - SX::mtimes(cAg(Slice(0, 6), Slice(6, int(Ag.cols()))), cvju))),
                        cvju)});

    SX cq0(model.nq);
    pinocchio::casadi::copy(q0_vec, cq0);

    Function L("L",
               {cx, cu},
               {1e-3 * SX::sumsqr(cvju) +
                1e-4 * SX::sumsqr(cf) +
                1e-4 * SX::sumsqr(ctau) +
                1e1 * SX::sumsqr(cqj - cq0(Slice(7, nq)))});

    Function Phi("Phi",
                 {cx},
                 {1e2 * SX::sumsqr(cqj - cq0(Slice(7, nq)))});

    Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 1;

    using namespace legged;
    using namespace constraints;

    shared_ptr<GeneralProblemData> gp_data = make_shared<GeneralProblemData>(Fint, Fdif, F, L, Phi);

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> friction_cone_constraint_builder =
        std::make_shared<FrictionConeConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> velocity_constraint_builder =
        std::make_shared<VelocityConstraintBuilder<LeggedRobotProblemData>>();

    std::shared_ptr<ConstraintBuilder<LeggedRobotProblemData>> contact_constraint_builder =
        std::make_shared<ContactConstraintBuilder<LeggedRobotProblemData>>();

    vector<shared_ptr<ConstraintBuilder<LeggedRobotProblemData>>> builders = {friction_cone_constraint_builder, velocity_constraint_builder, contact_constraint_builder};

    shared_ptr<LeggedRobotProblemData> legged_problem_data = make_shared<LeggedRobotProblemData>(gp_data, surfaces, contact_sequence, si, make_shared<Model>(model), make_shared<Data>(data), make_shared<ADModel>(cmodel), make_shared<ADData>(cdata), bot.getEndEffectors(), cx, cu, cdt, 20);

    std::vector<opt::ConstraintData> constraint_datas;
    for (auto builder : builders)
    {
        std::cout << "Building constraint" << std::endl;
        opt::ConstraintData some_data;
        builder->BuildConstraint(*legged_problem_data, 0, some_data);
        constraint_datas.push_back(some_data);
    }

    TrajectoryOpt<LeggedRobotProblemData, PseudospectralSegment> traj(legged_problem_data, builders, opts);

    DM X0 = DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    traj.init_finite_elements(1, X0);

    auto sol = traj.optimize();
    cout << sol << endl;

    return 0;
}