#include <pinocchio/autodiff/casadi.hpp>

using namespace casadi;

//quat = [qr qx qy qz]
SX quat_assign_if_negative(SX R, int i)
{
    auto quat = SX::zeros(4);
    auto j = (i + 1) % 3;
    auto k = (j + 1) % 3;
    auto t = sqrt(R(i,i) - R(j,j) - R(k,k) + 1.0);
    quat(i) = t / 2.0;
    t = 0.5 / t;
    quat(0) = (R(k,j) - R(j,k)) * t;
    quat(j) = (R(j,i) + R(i,j)) * t;
    quat(k) = (R(k,i) + R(i,k)) * t;
    return quat;
}

SX quat_assign_if_positive(SX R, SX t)
{
    auto quat = SX::zeros(4);
    auto t2 = sqrt(t + 1.0);
    quat(0) = t2 / 2.0;
    t2 = 0.5 / t2;
    quat(1) = (R(2,1) - R(1,2)) * t2;
    quat(2) = (R(0,2) - R(2,0)) * t2;
    quat(3) = (R(1,0) - R(0,1)) * t2;
    return quat;
}

SX lie_group_int(SX x, SX dx, SX dt)
{
    std::cout << "Running lie_group_int" << std::endl;
    auto pos = x(Slice(0, 3));
    auto quat = x(Slice(3, 7));
    auto r0c0 = pow(quat(0), 2) + pow(quat(1), 2) - pow(quat(2), 2) - pow(quat(3), 2);
    auto r0c1 = 2 * (quat(1) * quat(2) - quat(0) * quat(3));
    auto r0c2 = 2 * (quat(0) * quat(2) + quat(1) * quat(3));
    auto r1c0 = 2 * (quat(0) * quat(3) + quat(1) * quat(2));
    auto r1c1 = pow(quat(0), 2) - pow(quat(1), 2) + pow(quat(2), 2) - pow(quat(3), 2);
    auto r1c2 = 2 * (quat(2) * quat(3) - quat(0) * quat(1));
    auto r2c0 = 2 * (quat(1) * quat(3) - quat(0) * quat(2));
    auto r2c1 = 2 * (quat(0) * quat(1) + quat(2) * quat(3));
    auto r2c2 = pow(quat(0), 2) - pow(quat(1), 2) - pow(quat(2), 2) + pow(quat(3), 2);
    auto R = vertcat(horzcat(r0c0, r0c1, r0c2), horzcat(r1c0, r1c1, r1c2), horzcat(r2c0, r2c1, r2c2));
    auto M0 = vertcat(horzcat(R, pos), horzcat(SX(0), SX(0), SX(0), SX(1)));

    std::cout << "M0: " << M0.size() << std::endl;

    auto vel = dx(Slice(0, 3)) * dt;
    auto omega = dx(Slice(3, 6)) * dt;

    auto t2 = pow(omega(0), 2) + pow(omega(1), 2) + pow(omega(2), 2);
    auto t = sqrt(t2 + eps);
    auto ct = cos(t);
    auto st = sin(t);
    auto inv_t2 = 1 / t2;

    auto alpha_wxv = if_else(t < eps,
                             (1 / 2) - (t2 / 24),
                             (1 - ct) * inv_t2);

    auto alpha_v = if_else(t < eps,
                           1 - (t2 / 6),
                           st / t);

    auto alpha_w = if_else(t < eps,
                           (1 / 6) - (t2 / 120),
                           (1 - alpha_v) * inv_t2);

    auto diagonal_term = if_else(t < eps,
                                 1 - (t2 / 2),
                                 ct);

    auto trans = (alpha_v * vel + (alpha_w * dot(omega, vel) * omega + alpha_wxv * cross(omega, vel)));

    auto rot = alpha_wxv * mtimes(omega, omega.T());

    rot(0, 1) -= alpha_v * omega(2);
    rot(1, 0) += alpha_v * omega(2);
    rot(0, 2) += alpha_v * omega(1);
    rot(2, 0) -= alpha_v * omega(1);
    rot(1, 2) -= alpha_v * omega(0);
    rot(2, 1) += alpha_v * omega(0);
    
    rot(0, 0) += diagonal_term;
    rot(1, 1) += diagonal_term;
    rot(2, 2) += diagonal_term;

    auto exp6 = vertcat(horzcat(rot, trans), horzcat(SX(0), SX(0), SX(0), SX(1)));

    auto M1 = mtimes(M0, exp6);

    std::cout << "M1: " << M1.size() << std::endl;

    auto R1 = M1(Slice(0, 3), Slice(0, 3));

    auto pos_res = M1(Slice(0, 3), Slice(3, 4));

    auto tr = trace(R1);
    auto quat_res = if_else(tr > 0,
                            quat_assign_if_positive(R1, tr),
                            SX::vertcat({0.25, 0.25, 0.25, 0.25}));

    auto dot_product = dot(quat_res, quat);
    for (int i = 0; i < 4; ++i)
    {
        quat_res(i) = if_else(dot_product < 0, -quat_res(i), quat_res(i));
    }

    // auto N2 = pow(quat_res(0), 2) + pow(quat_res(1), 2) + pow(quat_res(2), 2) + pow(quat_res(3), 2);
    // auto alpha = (SX(3) - N2) / SX(2);
    // quat_res *= alpha;

    return vertcat(pos_res, quat_res);
}

int main()
{
    SX x = SX::sym("x", 7);
    SX dx = SX::sym("dx", 6);
    SX dt = SX::sym("dt");
    // Initialize x, dx, and dt
    DM x_val = DM::vertcat({DM::zeros(3), DM::vertcat({1, 0, 0, 0})}); // position at origin, quaternion for no rotation
    DM dx_val = DM::vertcat({DM::ones(3), DM::ones(3)});               // velocity and angular velocity of 1
    double dt_val = 1.0;                                               // time step of 1 second

    // Create a function for the lie group integrator
    Function lie_group_int_func = Function("lie_group_int", {x, dx, dt}, {lie_group_int(x, dx, dt)});
    DM result = lie_group_int_func(DMVector{x_val, dx_val, dt_val})[0];
    std::cout << "Result: " << result << std::endl;
    std::cout << "Norm Quaternion: " << DM::norm_2(result(Slice(3, 7))) << std::endl;

    return 0;
}
