#include "galileo/model/LeggedBody.h"
#include "galileo/variables/TrajectoryOpt.h"
#include <string>

#include <pinocchio/parsers/urdf.hpp>
#include <Eigen/Dense>

#include <pinocchio/autodiff/casadi.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>

#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/centroidal-derivatives.hpp>

using namespace galileo;

typedef double Scalar;
typedef casadi::SX ADScalar;

typedef pinocchio::ModelTpl<Scalar> Model;
typedef Model::Data Data;

typedef pinocchio::ModelTpl<ADScalar> ADModel;
typedef ADModel::Data ADData;

typedef Model::ConfigVectorType ConfigVector;
typedef Model::TangentVectorType TangentVector;

typedef ADModel::ConfigVectorType ConfigVectorAD;
typedef ADModel::TangentVectorType TangentVectorAD;

const std::string huron_location = "../resources/urdf/huron_cheat.urdf";

int main()
{
    double q0[] = {
        0, 0, 1.0627, 0, 0, 0, 1, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    // Map the array to Eigen matrix
    Eigen::Map<ConfigVector> q0_vec(q0, 19);

    std::vector<double> q0_stdvec = {
        0, 0, 1.0627, 0, 0, 0, 1, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365,
        0.0000, 0.0000, 0.0000, -0.3207, 0.7572, -0.4365, 0.0000};

    Model model;
    pinocchio::urdf::buildModel(huron_location, model);
    Data data(model);

    ADModel cmodel = model.cast<ADScalar>();
    ADData cdata(cmodel);

    pinocchio::computeTotalMass(model, data);
    pinocchio::framesForwardKinematics(model, data, q0_vec);

    auto mass = data.mass[0];
    auto g = casadi::SX::zeros(3, 1);
    g(2) = 9.81;
    auto nq = model.nq;
    auto nv = model.nv;
    std::shared_ptr<variables::States> si = std::make_shared<variables::States>(nq, nv);

    casadi::SX cx = casadi::SX::sym("x", si->nx);
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
    auto cf = si->get_f(cu);
    auto ctau = si->get_tau(cu);
    auto cvju = si->get_vju(cu);

    std::cout << "ch.size: " << ch.size() << std::endl;
    std::cout << "ch_d.size: " << ch_d.size() << std::endl;
    std::cout << "cdh.size: " << cdh.size() << std::endl;
    std::cout << "cdh_d.size: " << cdh_d.size() << std::endl;
    std::cout << "cq.size: " << cq.size() << std::endl;
    std::cout << "cq_d.size: " << cq_d.size() << std::endl;
    std::cout << "cqj.size: " << cqj.size() << std::endl;
    std::cout << "cv.size: " << cv.size() << std::endl;
    std::cout << "cv_d.size: " << cv_d.size() << std::endl;
    std::cout << "cvj.size: " << cvj.size() << std::endl;
    std::cout << "cf.size: " << cf.size() << std::endl;
    std::cout << "ctau.size: " << ctau.size() << std::endl;
    std::cout << "cvju.size: " << cvju.size() << std::endl;

    ConfigVectorAD cq_AD(model.nq);
    for (Eigen::DenseIndex k = 0; k < model.nq; ++k)
    {
        cq_AD[k] = cq(k);
    }

    TangentVectorAD cv_AD(model.nv);
    for (Eigen::DenseIndex k = 0; k < model.nv; ++k)
    {
        cv_AD[k] = cv(k);
    }

    TangentVectorAD cq_d_AD(model.nv);
    for (Eigen::DenseIndex k = 0; k < model.nv; ++k)
    {
        cq_d_AD[k] = cq_d(k);
    }

    pinocchio::centerOfMass(cmodel, cdata, cq_AD, false);
    pinocchio::computeCentroidalMap(cmodel, cdata, cq_AD);
    pinocchio::forwardKinematics(cmodel, cdata, cq_AD, cv_AD);
    pinocchio::updateFramePlacements(cmodel, cdata);

    auto intres = pinocchio::integrate(cmodel, cq_AD, cq_d_AD * cdt);
    casadi::SX tmp1 = casadi::SX::zeros(intres.rows(), 1);

    pinocchio::casadi::copy(intres, tmp1);
    auto Fint = casadi::Function("Fint",
                                 casadi::SXVector{cx, cdx, cdt},
                                 casadi::SXVector{vertcat(ch_d,
                                                          cdh_d,
                                                          tmp1,
                                                          cv_d)});

    auto Ag = cdata.Ag;
    casadi::SX tmp2 = casadi::SX(casadi::Sparsity::dense(Ag.rows(), Ag.cols()));
    pinocchio::casadi::copy(Ag, tmp2);

    auto F = casadi::Function("F",
                              casadi::SXVector{cx, cu},
                              casadi::SXVector{vertcat(cdh,
                                                       (cf - mass * g) / mass,
                                                       ctau / mass,
                                                       cv,
                                                       casadi::SX::mtimes(casadi::SX::inv(tmp2(casadi::Slice(0, 6), casadi::Slice(0, 6))), (mass * ch - casadi::SX::mtimes(tmp2(casadi::Slice(0, 6), casadi::Slice(6, int(Ag.cols()))), cvju))),
                                                       cvju)});

    casadi::SX tmp3 = casadi::SX(casadi::Sparsity::dense(q0_vec.rows(), 1));
    pinocchio::casadi::copy(q0_vec, tmp3);

    // auto L = casadi::Function("L",
    //                           casadi::SXVector{cx, cu},
    //                           casadi::SXVector{1e-3 * casadi::SX::sumsqr(cvju) +
    //                                            1e-4 * casadi::SX::sumsqr(cf) +
    //                                            1e-4 * casadi::SX::sumsqr(ctau) +
    //                                            1e1 * casadi::SX::sumsqr(cqj - tmp3(casadi::Slice(7, nq)))});

    // /*Dummy terminal cost*/
    // auto Phi = casadi::Function("Phi",
    //                             casadi::SXVector{cx},
    //                             casadi::SXVector{1e2 * casadi::SX::sumsqr(cqj - tmp3(casadi::Slice(7, nq)))});

    auto L = casadi::Function("L",
                              casadi::SXVector{cx, cu},
                              casadi::SXVector{0.});

    /*Dummy terminal cost*/
    auto Phi = casadi::Function("Phi",
                                casadi::SXVector{cx},
                                casadi::SXVector{0.});

    std::cout << "Fint outputs: " << vertcat(ch_d, cdh_d, tmp1, cv_d) << std::endl;

    std::cout << "F outputs" << vertcat(cdh, (cf - mass * g) / mass, ctau / mass, cv, casadi::SX::mtimes(casadi::SX::inv(tmp2(casadi::Slice(0, 6), casadi::Slice(0, 6))), (mass * ch - casadi::SX::mtimes(tmp2(casadi::Slice(0, 6), casadi::Slice(6, int(Ag.cols()))), cvju))), cvju) << std::endl;

    casadi::Dict opts;
    opts["ipopt.linear_solver"] = "ma97";
    std::shared_ptr<variables::ProblemData> problem = std::make_shared<variables::ProblemData>(Fint, F, L, Phi);
    variables::TrajectoryOpt traj(opts, si, problem);

    casadi::DM X0 = casadi::DM::zeros(si->nx, 1);
    int j = 0;
    for (int i = si->nh + si->ndh; i < si->nh + si->ndh + si->nq; ++i)
    {
        X0(i) = q0_stdvec[j];
        ++j;
    }

    std::cout << "X0: " << X0 << std::endl;

    traj.init_finite_elements(2, X0);

    auto sol = traj.optimize();
    std::cout << sol << std::endl;

    return 0;
}