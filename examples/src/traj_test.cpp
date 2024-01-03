#include "traj_test.h"

int main()
{
    double q0[] = {
        0, 0, 1.0627, 1, 0, 0, 0, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    Model model;
    pinocchio::urdf::buildModel(huron_location, model);
    Data data(model);

    pinocchio::computeTotalMass(model, data);
    pinocchio::framesForwardKinematics(model, data, q0_vec);

    ADModel cmodel = model.cast<ADScalar>();
    ADData cdata(cmodel);

    auto mass = data.mass[0];
    auto g = SX::zeros(3, 1);
    g(2) = 9.81;
    auto nq = model.nq;
    auto nv = model.nv;
    shared_ptr<opt::States> si = make_shared<opt::States>(nq, nv);

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

    // ConfigVector q0test(model.nq);
    // q0test[0] = 0;
    // q0test[1] = 0;
    // q0test[2] = 1.627;
    // q0test[3] = 0;
    // q0test[4] = 0;
    // q0test[5] = 0;
    // q0test[6] = 1;

    // TangentVector dqtest(model.nv);
    // dqtest[0] = 0;
    // dqtest[1] = 0;
    // dqtest[2] = 0;
    // dqtest[3] = 0;
    // dqtest[4] = 0;
    // dqtest[5] = 0;

    // auto q_res_test = pinocchio::integrate(model, q0test, dqtest);
    // cout << "cpp integrate: \n"
    //      << q_res_test << endl;
    // ConfigVectorAD pq0res(model.nq);
    // pq0res[0] = 0.06342760080100329;
    // pq0res[1] = 0.27635775159811415;
    // pq0res[2] = 1.9207733728707947;
    // pq0res[3] = -0.05961344771568522;
    // pq0res[4] = -0.624812847137685;
    // pq0res[5] = 0.29887279867696365;
    // pq0res[6] = 0.7188394765969299;
    // cout << "python integrate: \n"
    //      << pq0res << endl;

    Function Fint("Fint",
                  {cx, cdx, cdt},
                  {vertcat(ch_d,
                           cdh_d,
                            cq_result, // returns different expression than python version, NO idea why
                        //    custom_fint(cx, cdx, cdt),
                           cv_d)});

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
    // opts["ipopt.warm_start_init_point"] = "yes";
    // opts["ipopt.warm_start_same_structure"] = "yes";
    // opts["ipopt.print_timing_statistics"] = "yes";
    opts["ipopt.linear_solver"] = "ma97";
    opts["ipopt.ma97_order"] = "metis";
    opts["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts["ipopt.max_iter"] = 1;
    // opts["ipopt.print_level"] = 5;
    shared_ptr<opt::GeneralProblemData> problem = make_shared<opt::GeneralProblemData>(Fint, Fdif, F, L, Phi);
    opt::TrajectoryOpt traj(opts, si, problem);

    DM X0 = DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    cout << "X0: " << X0 << endl;

    traj.init_finite_elements(1, X0);

    auto sol = traj.optimize();
    cout << sol << endl;

    return 0;
}