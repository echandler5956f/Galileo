#pragma once

#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace math
    {
        /*-----------------------------------------------------
        A Collection of commonly used quaternion operations.

        Specialized for numerical robustness and use with Casadi.
        -----------------------------------------------------*/

        /**
         * @brief  Performs the Lie group integration of a state.
         *
         * This function takes the current pose, represented by the position (pos),
         * quaternion (quat), as well as a twist, composed of velocity (vel), and angular velocity (omega), and integrates it over a
         * given time step (dt) using the Lie group integration method. The resulting pose is returned
         * as a concatenated vector of the new position (pos_new) and quaternion (quat_new).
         *
         * @tparam Scalar Casadi type
         * @param qb The position and quaternion [pos, quat]
         * @param vb The velocity and angular velocity [vel, omega]
         * @param dt The time step
         * @return Scalar new state
         */
        template <typename Scalar>
        Scalar lie_group_int(Scalar qb, Scalar vb, Scalar dt);

        /**
         * @brief Numerically stable quaternion exponentiation.
         *
         * Given a quaternion, this function computes the exponential of the quaternion.
         * The exponential of a quaternion is defined as [exp(a) * v * sin(norm(v)) / norm(v), exp(a) * cos(norm(v))],
         * where 'a' is the scalar part of the quaternion and 'v' is the vector part of the quaternion. It is the result
         * of rotating a unit vector by the angle 'norm(v)' around the axis 'v' over a unit time interval.
         *
         * @tparam Scalar Casadi type
         * @param quat Quaternion [v, a]
         * @return Scalar exponential quaternion
         */
        template <typename Scalar>
        Scalar quat_exp(Scalar quat);

        /**
         * @brief Computes the Lie group difference between two configurations.
         *
         * This function calculates the Lie group difference between two configurations
         * represented by `qb1` and `qb2`. It returns the velocity and angular velocity
         * required to move from `qb1` to `qb2` over a time step `dt`.
         *
         * @tparam Scalar Casadi type
         * @param qb1 The first state [pos, quat]
         * @param qb2 The second state [pos, quat]
         * @param dt Time step
         * @return Scalar The difference between the two states
         */
        template <typename Scalar>
        Scalar lie_group_diff(Scalar qb1, Scalar qb2, Scalar dt);

        /**
         * @brief Numerically stable quaternion logarithm.
         *
         * This function takes a quaternion as input and computes the logarithm of the quaternion.
         * The logarithm of a quaternion is a vector that represents the rotation axis and the
         * angle of the quaternion. It is the tangent space vector required to rotate from the
         * identity quaternion to the input quaternion over a unit time interval.
         *
         * @tparam Scalar Casadi type
         * @param quat Quaternion [v, a]
         * @return Scalar Logarithm of the quaternion
         */
        template <typename Scalar>
        Scalar quat_log(Scalar quat);

        /**
         * @brief Applies a quaternion rotation to a 3D vector.
         *
         * @tparam Scalar Casadi type
         * @param quat Quaternion [v, a]
         * @param vec3 3D vector [x, y, z]
         * @return Scalar Rotated 3D vector
         */
        template <typename Scalar>
        Scalar apply_quat(Scalar quat, Scalar vec3);

        /**
         * @brief Composition of two quaternions.
         *
         * @tparam Scalar Casadi type
         * @param quat1 The first quaternion [v, a]
         * @param quat2 The second quaternion [v, a]
         * @return Scalar The composition of the two quaternions
         */
        template <typename Scalar>
        Scalar quat_mult(Scalar quat1, Scalar quat2);

        /**
         * @brief Computes the inverse of a quaternion.
         *
         * @tparam Scalar Casadi type
         * @param quat The quaternion [v, a]
         * @return Scalar The inverse of the quaternion
         */
        template <typename Scalar>
        Scalar quat_inv(Scalar quat);

        /**
         * @brief Computes the Rodrigues rotation matrix.
         *
         * @tparam Scalar Casadi type
         * @param omega The rotation vector [wx, wy, wz]
         * @return Scalar The Rodrigues rotation matrix
         */
        template <typename Scalar>
        Scalar rodrigues(Scalar omega);

        /**
         * @brief Computes the distance between two quaternions.
         *
         * @tparam Scalar Casadi type
         * @param quat1 The first quaternion [v, a]
         * @param quat2 The second quaternion [v, a]
         * @return Scalar The distance between the two quaternions
         */
        template <typename Scalar>
        Scalar quat_distance(Scalar quat1, Scalar quat2);

        /**
         * @brief Spherical linear interpolation between two quaternions.
         * 
         * Implements the fast SLERP algorithm from 
         * "A Fast and Accurate Estimate for Slerp, www.geometrictools.com/Documentation/FastAndAccurateSlerp.pdf"
         * 
         * This implementation relies soley on addition and subtraction, so there is no costly branching or trigonometric evaluations.
         * This is probably overkill.
         * 
         * @tparam Scalar Casadi type
         * @param quat1 The first quaternion [v, a]
         * @param quat2 The second quaternion [v, a]
         * @param t Interpolation parameter
         * @return Scalar The interpolated quaternion
         */
        template <typename Scalar>
        Scalar quat_slerp(Scalar quat1, Scalar quat2, Scalar t);
    }
}
