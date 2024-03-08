#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/autodiff/casadi.hpp>
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <random>

using namespace casadi;

casadi::SX customFint(casadi::SX x, casadi::SX dx, casadi::SX dt)
{
    casadi::SX x_12 = x(0);
    casadi::SX x_13 = x(1);
    casadi::SX x_14 = x(2);
    casadi::SX x_15 = x(3);
    casadi::SX x_16 = x(4);
    casadi::SX x_17 = x(5);
    casadi::SX x_18 = x(6);
    casadi::SX x_19 = x(7);
    casadi::SX x_20 = x(8);
    casadi::SX x_21 = x(9);
    casadi::SX x_22 = x(10);
    casadi::SX x_23 = x(11);
    casadi::SX x_24 = x(12);
    casadi::SX x_25 = x(13);
    casadi::SX x_26 = x(14);
    casadi::SX x_27 = x(15);
    casadi::SX x_28 = x(16);
    casadi::SX x_29 = x(17);
    casadi::SX x_30 = x(18);
    casadi::SX dx_12 = dx(0);
    casadi::SX dx_13 = dx(1);
    casadi::SX dx_14 = dx(2);
    casadi::SX dx_15 = dx(3);
    casadi::SX dx_16 = dx(4);
    casadi::SX dx_17 = dx(5);
    casadi::SX dx_18 = dx(6);
    casadi::SX dx_19 = dx(7);
    casadi::SX dx_20 = dx(8);
    casadi::SX dx_21 = dx(9);
    casadi::SX dx_22 = dx(10);
    casadi::SX dx_23 = dx(11);
    casadi::SX dx_24 = dx(12);
    casadi::SX dx_25 = dx(13);
    casadi::SX dx_26 = dx(14);
    casadi::SX dx_27 = dx(15);
    casadi::SX dx_28 = dx(16);
    casadi::SX dx_29 = dx(17);

    casadi::SX c1 = (dx_12 * dt);
    casadi::SX c2 = (dx_15 * dt);
    casadi::SX c3 = (dx_16 * dt);
    casadi::SX c4 = (dx_17 * dt);
    casadi::SX c5 = 4.93038e-32;
    casadi::SX c6 = ((sq(c2) + (sq(c3) + sq(c4))) + c5);
    casadi::SX c7 = sqrt(c6);
    casadi::SX c8 = 0.00012207;
    casadi::SX c9 = (c7 < c8);
    casadi::SX c10 = 0.5;
    casadi::SX c11 = 24;
    casadi::SX c12 = 1;
    casadi::SX c13 = (if_else(c9, (c10 - (c6 / c11)), 0) + if_else((!c9), ((c12 - cos(c7)) / c6), 0));
    casadi::SX c14 = (dx_14 * dt);
    casadi::SX c15 = (dx_13 * dt);
    casadi::SX c16 = (c7 < c8);
    casadi::SX c17 = 120;
    casadi::SX c18 = (if_else(c16, (0.166667 - (c6 / c17)), 0) + if_else((!c16), (((c7 - sin(c7)) / c6) / c7), 0));
    casadi::SX c19 = ((c2 * c15) - (c3 * c1));
    casadi::SX c20 = ((c4 * c1) - (c2 * c14));
    casadi::SX c21 = ((c1 + (c13 * ((c3 * c14) - (c4 * c15)))) + (c18 * ((c3 * c19) - (c4 * c20))));
    casadi::SX c22 = ((c3 * c14) - (c4 * c15));
    casadi::SX c23 = ((c14 + (c13 * ((c2 * c15) - (c3 * c1)))) + (c18 * ((c2 * c20) - (c3 * c22))));
    casadi::SX c24 = ((c15 + (c13 * ((c4 * c1) - (c2 * c14)))) + (c18 * ((c4 * c22) - (c2 * c19))));
    casadi::SX c25 = ((x_16 * c23) - (x_17 * c24));
    casadi::SX c26 = (c25 + c25);
    casadi::SX c27 = ((x_15 * c24) - (x_16 * c21));
    casadi::SX c28 = (c27 + c27);
    casadi::SX c29 = ((x_17 * c21) - (x_15 * c23));
    casadi::SX c30 = (c29 + c29);
    casadi::SX c31 = (sq(c2) + (sq(c3) + sq(c4)));
    casadi::SX c32 = (c8 < c31);
    casadi::SX c33 = sqrt((c31 + c5));
    casadi::SX c34 = (c10 * c33);
    casadi::SX c35 = sin(c34);
    casadi::SX c36 = (c31 / 4);
    casadi::SX c37 = (c10 * ((c12 - (c36 / 6)) + (sq(c36) / c17)));
    casadi::SX c38 = (if_else(c32, (c35 * (c2 / c33)), 0) + if_else((!c32), (c37 * c2), 0));
    casadi::SX c39 = (c8 < c31);
    casadi::SX c40 = 2;
    casadi::SX c41 = (if_else(c39, cos(c34), 0) + if_else((!c39), ((c12 - (c36 / c40)) + (sq(c36) / c11)), 0));
    casadi::SX c42 = (c8 < c31);
    casadi::SX c43 = (if_else(c42, (c35 * (c4 / c33)), 0) + if_else((!c42), (c37 * c4), 0));
    casadi::SX c44 = (c8 < c31);
    casadi::SX c45 = (if_else(c44, (c35 * (c3 / c33)), 0) + if_else((!c44), (c37 * c3), 0));
    casadi::SX c46 = ((((x_18 * c38) + (x_15 * c41)) + (x_16 * c43)) - (x_17 * c45));
    casadi::SX c47 = ((((x_18 * c45) + (x_16 * c41)) + (x_17 * c38)) - (x_15 * c43));
    casadi::SX c48 = ((((x_18 * c43) + (x_17 * c41)) + (x_15 * c45)) - (x_16 * c38));
    casadi::SX c49 = ((((x_18 * c41) - (x_15 * c38)) - (x_16 * c45)) - (x_17 * c43));
    casadi::SX c50 = (((c46 * x_15) + (c47 * x_16)) + ((c48 * x_17) + (c49 * x_18)));
    casadi::SX c51 = 0;
    casadi::SX c52 = (c50 < c51);
    casadi::SX c53 = (if_else(c52, (-c46), 0) + if_else((!c52), c46, 0));
    casadi::SX c54 = (c50 < c51);
    casadi::SX c55 = (if_else(c54, (-c47), 0) + if_else((!c54), c47, 0));
    casadi::SX c56 = (c50 < c51);
    casadi::SX c57 = (if_else(c56, (-c48), 0) + if_else((!c56), c48, 0));
    casadi::SX c58 = (c50 < c51);
    casadi::SX c59 = (if_else(c58, (-c49), 0) + if_else((!c58), c49, 0));
    casadi::SX c60 = ((3 - ((sq(c53) + sq(c55)) + (sq(c57) + sq(c59)))) / c40);

    return casadi::SX::vertcat(
        {(((c21 + (x_18 * c26)) + ((x_16 * c28) - (x_17 * c30))) + x_12),
         (((c24 + (x_18 * c30)) + ((x_17 * c26) - (x_15 * c28))) + x_13),
         (((c23 + (x_18 * c28)) + ((x_15 * c30) - (x_16 * c26))) + x_14),
         (c53 * c60),
         (c55 * c60),
         (c57 * c60),
         (c59 * c60),
         (x_19 + (dx_18 * dt)),
         (x_20 + (dx_19 * dt)),
         (x_21 + (dx_20 * dt)),
         (x_22 + (dx_21 * dt)),
         (x_23 + (dx_22 * dt)),
         (x_24 + (dx_23 * dt)),
         (x_25 + (dx_24 * dt)),
         (x_26 + (dx_25 * dt)),
         (x_27 + (dx_26 * dt)),
         (x_28 + (dx_27 * dt)),
         (x_29 + (dx_28 * dt)),
         (x_30 + (dx_29 * dt))});
}

// quat = [qr qx qy qz]
SX quat_assign_if_negative(SX R, int i)
{
    SX quat = SX::zeros(4, 1);
    int j = (i + 1) % 3;
    int k = (j + 1) % 3;
    SX t = sqrt(R(i, i) - R(j, j) - R(k, k) + 1.0);
    quat(i) = t / 2.0;
    t = 0.5 / t;
    quat(0) = (R(k, j) - R(j, k)) * t;
    quat(j) = (R(j, i) + R(i, j)) * t;
    quat(k) = (R(k, i) + R(i, k)) * t;
    return quat;
}

SX quat_assign_if_positive(SX R, SX t)
{
    SX quat = SX::zeros(4, 1);
    SX t2 = sqrt(t + 1.0);
    quat(0) = t2 / 2.0;
    t2 = 0.5 / t2;
    quat(1) = (R(2, 1) - R(1, 2)) * t2;
    quat(2) = (R(0, 2) - R(2, 0)) * t2;
    quat(3) = (R(1, 0) - R(0, 1)) * t2;
    return quat;
}

unsigned long long factorial(int n)
{
    unsigned long long fact = 1;
    for (int i = 1; i <= n; ++i)
    {
        fact *= i;
    }
    return fact;
}

template <typename T>
void quat_to_rot(T quat, T &R_quat)
{

    auto a = quat(3);
    auto b = quat(0);
    auto c = quat(1);
    auto d = quat(2);

    auto s = 2 / (a * a + b * b + c * c + d * d);
    auto bs = b * s;
    auto cs = c * s;
    auto ds = d * s;
    auto ab = a * bs;
    auto ac = a * cs;
    auto ad = a * ds;
    auto bb = b * bs;
    auto bc = b * cs;
    auto bd = b * ds;
    auto cc = c * cs;
    auto cd = c * ds;
    auto dd = d * ds;

    R_quat(0, 0) = 1 - cc - dd;
    R_quat(0, 1) = bc - ad;
    R_quat(0, 2) = bd + ac;
    R_quat(1, 0) = bc + ad;
    R_quat(1, 1) = 1 - bb - dd;
    R_quat(1, 2) = cd - ab;
    R_quat(2, 0) = bd - ac;
    R_quat(2, 1) = cd + ab;
    R_quat(2, 2) = 1 - bb - cc;
}

void quat_to_rot_eigen(Eigen::VectorXd quat, Eigen::Matrix3d &R_quat)
{
    auto a = quat(3);
    auto b = quat(0);
    auto c = quat(1);
    auto d = quat(2);

    auto s = 2 / (a * a + b * b + c * c + d * d);
    auto bs = b * s;
    auto cs = c * s;
    auto ds = d * s;
    auto ab = a * bs;
    auto ac = a * cs;
    auto ad = a * ds;
    auto bb = b * bs;
    auto bc = b * cs;
    auto bd = b * ds;
    auto cc = c * cs;
    auto cd = c * ds;
    auto dd = d * ds;

    R_quat(0, 0) = 1 - cc - dd;
    R_quat(0, 1) = bc - ad;
    R_quat(0, 2) = bd + ac;
    R_quat(1, 0) = bc + ad;
    R_quat(1, 1) = 1 - bb - dd;
    R_quat(1, 2) = cd - ab;
    R_quat(2, 0) = bd - ac;
    R_quat(2, 1) = cd + ab;
    R_quat(2, 2) = 1 - bb - cc;
}

SX SE3_from_Rt(SX R, SX t)
{
    SX SE3 = SX::zeros(4, 4);
    SE3(Slice(0, 3), Slice(0, 3)) = R;
    SE3(Slice(0, 3), Slice(3, 4)) = t;
    SE3(3, 3) = 1;
    return SE3;
}

SX quat_norm(SX quat)
{
    return quat / sqrt(quat(0) * quat(0) + quat(1) * quat(1) + quat(2) * quat(2) + quat(3) * quat(3) + eps);
}

SX apply_quat(SX quat, SX vec3)
{
    SX imag = quat(Slice(0, 3));
    SX real = quat(3);
    return vec3 + 2 * cross(imag, cross(imag, vec3) + real * vec3);
}

SX omega_to_unit_quat(SX omega)
{
    // SX quat = SX::zeros(4, 1);
    // SX t = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2) + eps);
    // SX s = sin(t / 2) / t;
    // return SX::vertcat({s * omega, cos(t / 2)});
    return SX::vertcat({omega, 0});
}

SX quat_mult(SX quat1, SX quat2)
{
    SX real_1 = quat1(3);
    SX imag_1 = quat1(Slice(0, 3));

    SX real_2 = quat2(3);
    SX imag_2 = quat2(Slice(0, 3));

    SX real_res = real_1 * real_2 - dot(imag_1, imag_2);
    SX imag_res = real_1 * imag_2 + real_2 * imag_1 + cross(imag_1, imag_2);
    return SX::vertcat({imag_res, real_res});
}

SX quat_exp(SX quat)
{
    SX v = quat(Slice(0, 3));
    SX v_norm = sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2) + eps);
    SX a = quat(3);

    SX exp_a = exp(a);
    SX cos_norm_v = cos(v_norm);
    SX sin_norm_v = sin(v_norm);

    SX quat_new = SX::vertcat({exp_a * v * sin_norm_v / v_norm, exp_a * cos_norm_v});

    return quat_new / sqrt(quat_new(0) * quat_new(0) + quat_new(1) * quat_new(1) + quat_new(2) * quat_new(2) + quat_new(3) * quat_new(3) + eps);
}

SX rot_exp_taylor_series(SX rot, int degree)
{
    SX exp_a = SX::zeros(3);
    for (int i = 0; i < degree; i++)
    {
        exp_a = exp_a + (mpower(rot, i) / factorial(i));
    }
    return exp_a;
}

SX rodrigues(SX omega)
{
    SX theta = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2) + eps);
    return SX::eye(3) + ((1 - cos(theta)) / (theta * theta)) * skew(omega) + ((theta - sin(theta)) / (theta * theta * theta)) * mpower(skew(omega), 2);
}

SX lie_group_int(SX x, SX dx, SX dt)
{
    SX pos = x(Slice(0, 3));
    SX quat = x(Slice(3, 7));

    SX vel = dx(Slice(0, 3)) * dt;
    SX omega = dx(Slice(3, 6)) * dt;
    SX omega_quat = omega_to_unit_quat(omega / 2);
    SX exp_omega_quat = quat_exp(omega_quat);
    SX quat_new = quat_norm(quat_mult(quat, exp_omega_quat));
    SX pos_new = pos + apply_quat(quat, SX::mtimes(rodrigues(omega), vel));
    SX joint_pos_new = x(Slice(7, x.size1())) + dx(Slice(6, dx.size1())) * dt;

    return vertcat(pos_new, quat_new, joint_pos_new);
}

SX alt_int(SX x, SX dx, SX dt, SX &R_new)
{
    SX pos = x(Slice(0, 3));
    SX quat = x(Slice(3, 7));

    SX vel = dx(Slice(0, 3)) * dt;
    SX omega = dx(Slice(3, 6)) * dt;

    SX rot_quat = SX::eye(3);
    quat_to_rot(quat, rot_quat);
    SX exp_omega_rot = rot_exp_taylor_series(skew(omega), 8);
    R_new = casadi::SX::mtimes(rot_quat, exp_omega_rot);
    SX pos_new = pos + apply_quat(quat, SX::mtimes(rodrigues(omega), vel));

    return pos_new;
}

int main()
{

    auto model = pinocchio::Model();
    pinocchio::urdf::buildModel("../../resources/go1/urdf/go1.urdf", pinocchio::JointModelFreeFlyer(), model);
    auto data = pinocchio::Data(model);

    SX cx = SX::sym("x", model.nq, 1);
    SX cdx = SX::sym("dx", model.nv, 1);
    SX cdt = SX::sym("dt", 1, 1);

    Function baseline_func = Function("Fint",
                                      {cx, cdx, cdt},
                                      {customFint(cx, cdx, cdt)});

    Function baseline_jac_func = Function("baseline_jac",
                                          {cx, cdx, cdt},
                                          {jacobian(customFint(cx, cdx, cdt), vertcat(SXVector{cx, cdx, cdt}))});

    Function lie_group_int_func = Function("lie_group_int",
                                           {cx, cdx, cdt},
                                           {lie_group_int(cx, cdx, cdt)});

    Function lie_group_int_jac_func = Function("lie_group_int_jac",
                                               {cx, cdx, cdt},
                                               {jacobian(lie_group_int(cx, cdx, cdt), vertcat(SXVector{cx, cdx, cdt}))});

    std::cout << "Comparing the two functions..." << std::endl;

    /*Compare the two functions*/
    int num_samples = 10000;
    DM error = DM::zeros(model.nq, 1);
    DM jac_error = DM::zeros(model.nq, model.nq + model.nv + 1);
    std::uniform_real_distribution<double> unif_pos(-1.0, 1.0);
    std::uniform_real_distribution<double> unif_time(0.0, 1.0);
    std::default_random_engine re;

    for (int i = 0; i < num_samples; ++i)
    {
        Eigen::VectorXd x = pinocchio::randomConfiguration(model);
        x[0] = unif_pos(re);
        x[1] = unif_pos(re);
        x[2] = unif_pos(re);
        Eigen::VectorXd dx = Eigen::VectorXd::Random(model.nv);
        DM dt = unif_time(re);
        DM x_dm(x.size(), 1);
        DM dx_dm(dx.size(), 1);
        pinocchio::casadi::copy(x, x_dm);
        pinocchio::casadi::copy(dx, dx_dm);

        DM baseline_res = baseline_func(DMVector{x_dm, dx_dm, dt}).at(0);
        DM lie_group_int_res = lie_group_int_func(DMVector{x_dm, dx_dm, dt}).at(0);
        DM baseline_jac_res = baseline_jac_func(DMVector{x_dm, dx_dm, dt}).at(0);
        DM lie_group_int_jac_res = lie_group_int_jac_func(DMVector{x_dm, dx_dm, dt}).at(0);

        // std::cout << "dx: " << dx_dm(Slice(0, 6)) << std::endl;
        // std::cout << "Baseline: " << baseline_res << std::endl;
        // std::cout << "Lie Group Int: " << lie_group_int_res << std::endl;
        // std::cout << "Baseline Jacobian: " << baseline_jac_res << std::endl;
        // std::cout << "Lie Group Int Jacobian: " << lie_group_int_jac_res << std::endl;
        // std::cout << std::endl;

        error += pow(baseline_res - lie_group_int_res, 2);
        jac_error += pow(baseline_jac_res - lie_group_int_jac_res, 2);
    }
    error = sqrt(error / num_samples);
    std::cout << "Error: " << error << std::endl;
    jac_error = sqrt(jac_error / num_samples);
    std::cout << "Jacobian Error: " << jac_error << std::endl;

    return 0;
}
