#include "galileo/math/LieAlgebra.h"

namespace galileo
{
    namespace math
    {
        template <typename Scalar>
        Scalar lie_group_int(Scalar qb, Scalar vb, Scalar dt)
        {
            assert(qb.size1() == 7);
            assert(vb.size1() == 6);
            assert(qb.size2() == 1);
            assert(vb.size2() == 1);

            Scalar pos = qb(casadi::Slice(0, 3));
            Scalar quat = qb(casadi::Slice(3, 7));

            Scalar vel = vb(casadi::Slice(0, 3)) * dt;
            Scalar omega = vb(casadi::Slice(3, 6)) * dt;
            Scalar omega_quat = Scalar::vertcat({omega / 2, 0});
            Scalar exp_omega_quat = quat_exp(omega_quat);
            Scalar quat_new = quat_mult(quat, exp_omega_quat);
            quat_new = quat_new / sqrt(quat_new(0) * quat_new(0) + quat_new(1) * quat_new(1) + quat_new(2) * quat_new(2) + quat_new(3) * quat_new(3) + casadi::eps * casadi::eps);
            Scalar pos_new = pos + apply_quat(quat, Scalar::mtimes(rodrigues(omega), vel));
            return vertcat(pos_new, quat_new);
        }

        template casadi::DM lie_group_int<casadi::DM>(casadi::DM qb, casadi::DM vb, casadi::DM dt);
        template casadi::SX lie_group_int<casadi::SX>(casadi::SX qb, casadi::SX vb, casadi::SX dt);
        template casadi::MX lie_group_int<casadi::MX>(casadi::MX qb, casadi::MX vb, casadi::MX dt);

        template <typename Scalar>
        Scalar quat_exp(Scalar quat)
        {
            assert(quat.size1() == 4);
            assert(quat.size2() == 1);

            Scalar v = quat(casadi::Slice(0, 3));
            Scalar v_norm = sqrt(v(0) * v(0) + v(1) * v(1) + v(2) * v(2) + casadi::eps * casadi::eps);
            Scalar a = quat(3);

            Scalar exp_a = exp(a);
            Scalar cos_norm_v = cos(v_norm);
            Scalar sin_norm_v = sin(v_norm);

            Scalar quat_new = Scalar::vertcat({exp_a * v * sin_norm_v / v_norm, exp_a * cos_norm_v});

            return quat_new / sqrt(quat_new(0) * quat_new(0) + quat_new(1) * quat_new(1) + quat_new(2) * quat_new(2) + quat_new(3) * quat_new(3) + casadi::eps * casadi::eps);
        }

        template casadi::DM quat_exp<casadi::DM>(casadi::DM quat);
        template casadi::SX quat_exp<casadi::SX>(casadi::SX quat);
        template casadi::MX quat_exp<casadi::MX>(casadi::MX quat);

        template <typename Scalar>
        Scalar lie_group_diff(Scalar qb1, Scalar qb2, Scalar dt)
        {
            assert(qb1.size1() == 7);
            assert(qb2.size1() == 7);
            assert(qb1.size2() == 1);
            assert(qb2.size2() == 1);

            Scalar pos1 = qb1(casadi::Slice(0, 3));
            Scalar quat1 = qb1(casadi::Slice(3, 7));
            Scalar pos2 = qb2(casadi::Slice(0, 3));
            Scalar quat2 = qb2(casadi::Slice(3, 7));

            Scalar quat_diff = quat_mult(quat_inv(quat1), quat2);
            Scalar omega_quat = quat_log(quat_diff);
            Scalar omega = omega_quat(casadi::Slice(0, 3));

            Scalar pos_diff = pos2 - pos1;
            Scalar vel = Scalar::mtimes(inv(rodrigues(omega)), apply_quat(quat_inv(quat1), pos_diff));

            return Scalar::vertcat({vel / dt, omega / dt});
        }

        template casadi::DM lie_group_diff<casadi::DM>(casadi::DM qb1, casadi::DM qb2, casadi::DM dt);
        template casadi::SX lie_group_diff<casadi::SX>(casadi::SX qb1, casadi::SX qb2, casadi::SX dt);
        template casadi::MX lie_group_diff<casadi::MX>(casadi::MX qb1, casadi::MX qb2, casadi::MX dt);

        template <typename Scalar>
        Scalar quat_log(Scalar quat)
        {
            assert(quat.size1() == 4);
            assert(quat.size2() == 1);

            Scalar squared_n = quat(0) * quat(0) + quat(1) * quat(1) + quat(2) * quat(2);
            Scalar n = sqrt(squared_n);
            Scalar w = quat(3);

            return Scalar::if_else(n < casadi::eps,
                                   ((2. / w) - (2. / 3.) * (squared_n / (w * w * w))) * quat(casadi::Slice(0, 3)),
                                   4. * atan(n / (w + sqrt(w * w + squared_n + casadi::eps * casadi::eps))) * quat(casadi::Slice(0, 3)) / n);
        }

        template casadi::DM quat_log<casadi::DM>(casadi::DM quat);
        template casadi::SX quat_log<casadi::SX>(casadi::SX quat);
        template casadi::MX quat_log<casadi::MX>(casadi::MX quat);

        template <typename Scalar>
        Scalar apply_quat(Scalar quat, Scalar vec3)
        {
            assert(quat.size1() == 4);
            assert(quat.size2() == 1);
            assert(vec3.size1() == 3);
            assert(vec3.size2() == 1);

            Scalar imag = quat(casadi::Slice(0, 3));
            Scalar real = quat(3);
            return vec3 + 2 * cross(imag, cross(imag, vec3) + real * vec3);
        }

        template casadi::DM apply_quat<casadi::DM>(casadi::DM quat, casadi::DM vec3);
        template casadi::SX apply_quat<casadi::SX>(casadi::SX quat, casadi::SX vec3);
        template casadi::MX apply_quat<casadi::MX>(casadi::MX quat, casadi::MX vec3);

        template <typename Scalar>
        Scalar quat_mult(Scalar quat1, Scalar quat2)
        {
            assert(quat1.size1() == 4);
            assert(quat1.size2() == 1);
            assert(quat2.size1() == 4);
            assert(quat2.size2() == 1);

            Scalar real_1 = quat1(3);
            Scalar imag_1 = quat1(casadi::Slice(0, 3));

            Scalar real_2 = quat2(3);
            Scalar imag_2 = quat2(casadi::Slice(0, 3));

            Scalar real_res = real_1 * real_2 - dot(imag_1, imag_2);
            Scalar imag_res = real_1 * imag_2 + real_2 * imag_1 + cross(imag_1, imag_2);
            return Scalar::vertcat({imag_res, real_res});
        }

        template casadi::DM quat_mult<casadi::DM>(casadi::DM quat1, casadi::DM quat2);
        template casadi::SX quat_mult<casadi::SX>(casadi::SX quat1, casadi::SX quat2);
        template casadi::MX quat_mult<casadi::MX>(casadi::MX quat1, casadi::MX quat2);

        template <typename Scalar>
        Scalar quat_inv(Scalar quat)
        {
            assert(quat.size1() == 4);
            assert(quat.size2() == 1);

            Scalar quat_norm = sqrt(quat(0) * quat(0) + quat(1) * quat(1) + quat(2) * quat(2) + quat(3) * quat(3) + casadi::eps * casadi::eps);
            return Scalar::vertcat({-quat(casadi::Slice(0, 3)) / quat_norm, quat(3) / quat_norm});
        }

        template casadi::DM quat_inv<casadi::DM>(casadi::DM quat);
        template casadi::SX quat_inv<casadi::SX>(casadi::SX quat);
        template casadi::MX quat_inv<casadi::MX>(casadi::MX quat);

        template <typename Scalar>
        Scalar rodrigues(Scalar omega)
        {
            assert(omega.size1() == 3);
            assert(omega.size2() == 1);

            Scalar theta = sqrt(omega(0) * omega(0) + omega(1) * omega(1) + omega(2) * omega(2) + casadi::eps * casadi::eps);
            return Scalar::eye(3) + ((1 - cos(theta)) / (theta * theta)) * skew(omega) + ((theta - sin(theta)) / (theta * theta * theta)) * mpower(skew(omega), 2);
        }

        template casadi::DM rodrigues<casadi::DM>(casadi::DM omega);
        template casadi::SX rodrigues<casadi::SX>(casadi::SX omega);
        template casadi::MX rodrigues<casadi::MX>(casadi::MX omega);

        template <typename Scalar>
        Scalar quat_distance(Scalar quat1, Scalar quat2)
        {
            assert(quat1.size1() == 4);
            assert(quat1.size2() == 1);
            assert(quat2.size1() == 4);
            assert(quat2.size2() == 1);

            Scalar real_1 = quat1(3);
            Scalar imag_1 = quat1(casadi::Slice(0, 3));

            Scalar real_2 = quat2(3);
            Scalar imag_2 = quat2(casadi::Slice(0, 3));

            return real_1 * imag_2 - real_2 * imag_1 + cross(imag_1, imag_2);
        }

        template casadi::DM quat_distance<casadi::DM>(casadi::DM quat1, casadi::DM quat2);
        template casadi::SX quat_distance<casadi::SX>(casadi::SX quat1, casadi::SX quat2);
        template casadi::MX quat_distance<casadi::MX>(casadi::MX quat1, casadi::MX quat2);
    }
}