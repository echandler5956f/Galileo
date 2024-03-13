#pragma once

#include "galileo/legged-model/ContactSequence.h"
#include <pinocchio/algorithm/center-of-mass.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <Eigen/LU>
#include <bitset>

namespace galileo
{
    namespace legged
    {
        /**
         * @brief A class for holding the robot model.
         *
         */
        class LeggedBody
        {
        public:
            /**
             * @brief Construct a new Legged Body object.
             *
             * @param location The location of the URDF file.
             * @param num_ees The number of end effectors.
             * @param end_effector_names The string IDs that correspond to the pinocchio end effector frames.
             *
             */
            LeggedBody(const std::string location, const int num_ees, const std::string end_effector_names[]);

            /**
             * @brief Provide the string IDs that correspond to the pinocchio end effector frames.
             *
             * @param ee_names The string IDs that correspond to the pinocchio end effector frames.
             */
            void setEndEffectors(const std::vector<std::string> &ee_names);

            /**
             * @brief Generate combinations of contacts.
             *
             */
            void generateContactCombination();

            /**
             * @brief Create the generalized dynamics, fint, and fdiff functions.
             *
             */
            void createGeneralFunctions();

            /**
             * @brief Create the phase-invariant centroidal momentum dynamics (summed wrenches).
             *
             */
            void createGeneralDynamics();

            /**
             * @brief Create the state integrator function to translate from the tangent space to the state space.
             *
             */
            void createFint();

            /**
             * @brief Performs the Lie group integration for the LeggedBody class.
             *
             * This function takes the current pose of the LeggedBody, represented by the position (pos),
             * quaternion (quat), velocity (vel), and angular velocity (omega), and integrates it over a
             * given time step (dt) using the Lie group integration method. The resulting pose is returned
             * as a concatenated vector of the new position (pos_new) and quaternion (quat_new).
             *
             * @param qb The current pose vector of the LeggedBody, consisting of the position (pos) and quaternion (quat).
             * @param vb The current velocity vector of the LeggedBody, consisting of the velocity (vel) and angular velocity (omega).
             * @param dt The time step for integration.
             * @return The new pose vector of the LeggedBody, consisting of the new position (pos_new) and quaternion (quat_new).
             */
            casadi::SX lie_group_int(casadi::SX qb, casadi::SX vb, casadi::SX dt);

            /**
             * @brief Numerically stable quaternion exponentiation
             *
             * Given a quaternion, this function computes the exponential of the quaternion.
             * The exponential of a quaternion is defined as [exp(a) * v * sin(norm(v)) / norm(v), exp(a) * cos(norm(v))],
             * where 'a' is the scalar part of the quaternion and 'v' is the vector part of the quaternion. It is the result
             * of rotating a unit vector by the angle 'norm(v)' around the axis 'v' over a unit time interval.
             *
             * @param quat The input quaternion.
             * @return The exponential of the input quaternion.
             */
            casadi::SX quat_exp(casadi::SX quat);

            /**
             * @brief Create the state derivative function to translate from the state space to the tangent space.
             *
             */
            void createFdiff();

            /**
             * @brief Computes the Lie group difference between two configurations.
             *
             * This function calculates the Lie group difference between two configurations
             * represented by `qb1` and `qb2`. It returns the velocity and angular velocity
             * of the difference in the form of a `casadi::SX` vector.
             *
             * @param qb1 The first configuration vector.
             * @param qb2 The second configuration vector.
             * @param dt The time step.
             * @return The Lie group difference between the two configurations.
             */
            casadi::SX lie_group_diff(casadi::SX qb1, casadi::SX qb2, casadi::SX dt);

            /**
             * @brief Numerically stable quaternion logarithm
             *
             * This function takes a quaternion as input and computes the logarithm of the quaternion.
             * The logarithm of a quaternion is a vector that represents the rotation axis and the
             * angle of the quaternion. It is the tangent space vector required to rotate from the
             * identity quaternion to the input quaternion over a unit time interval.
             *
             * @param quat The input quaternion.
             * @return The logarithm of the quaternion.
             */
            casadi::SX quat_log(casadi::SX quat);

            /**
             * @brief Apply a quaternion to a 3D vector.
             *
             * @param quat Orientation quaternion to apply.
             * @param vec3 Vector to rotate.
             * @return casadi::SX The rotated vector.
             */
            casadi::SX apply_quat(casadi::SX quat, casadi::SX vec3);

            /**
             * @brief Multiply two quaternions.
             *
             * @param quat1 The first quaternion.
             * @param quat2 The second quaternion.
             * @return casadi::SX The product of the two quaternions.
             */
            casadi::SX quat_mult(casadi::SX quat1, casadi::SX quat2);

            /**
             * @brief Invert a quaternion.
             *
             * @param quat The quaternion to invert.
             * @return casadi::SX The inverted quaternion.
             */
            casadi::SX quat_inv(casadi::SX quat);

            /**
             * @brief Calculates the Rodrigues rotation matrix for a given angular velocity vector.
             *
             * @param omega The angular velocity vector.
             * @return The Rodrigues rotation matrix.
             */
            casadi::SX rodrigues(casadi::SX omega);

            /**
             * @brief Create the state error function.
             *
             */
            void createErrorFunction();

            /**
             * @brief Calculate the distance between two quaternions.
             *
             * The resulting 3x1 vector is a measure of how much one quaternion is rotated from another.
             * If this vector is zero, it indicates that the measured and desired frames are aligned,
             * meaning the orientations are the same. If not, the vector gives the axis and magnitude
             * of rotation needed to align the measured quaternion with the desired one.
             *
             * @param quat1 The first quaternion.
             * @param quat2 The second quaternion.
             *
             * @return casadi::SX The distance between the two quaternions.
             */
            casadi::SX quat_distance(casadi::SX quat1, casadi::SX quat2);

            /**
             * @brief Create the mode dynamics for each mode using the General Dynamics.
             *
             *  @param print_ees_info Whether to print the end effector contact information.
             */
            void fillModeDynamics(bool print_ees_info = false);

            /**
             * @brief Get the Contact Combination object from a binary combination mask.
             *
             * @param contact_mask The binary combination mask.
             * @return contact::ContactCombination The contact combination.
             */
            contact::ContactCombination getContactCombination(int contact_mask);

            /**
             * @brief Get the End Effectors.
             *
             */
            contact::RobotEndEffectors getEndEffectors();

            /**
             * @brief The pinocchio model of the robot.
             *
             */
            opt::Model model;

            /**
             * @brief The pinocchio data of the robot.
             *
             */
            opt::Data data;

            /**
             * @brief The symbolic model of the robot.
             *
             */
            opt::ADModel cmodel;

            /**
             * @brief The symbolic data of the robot.
             *
             */
            opt::ADData cdata;

            /**
             * @brief The state indices helper for the robot.
             *
             */
            std::shared_ptr<opt::LeggedRobotStates> si;

            /**
             * @brief The contact sequence of the robot.
             *
             */
            std::shared_ptr<contact::ContactSequence> contact_sequence;

            /**
             * @brief The general dynamics function.
             *
             */
            casadi::Function general_dynamics;

            /**
             * @brief The state integrator function.
             *
             */
            casadi::Function fint;

            /**
             * @brief The state derivative function.
             *
             */
            casadi::Function fdiff;

            /**
             * @brief The state error function
             *
             */
            casadi::Function f_state_error;

            /**
             * @brief The symbolic state.
             *
             */
            casadi::SX cx;

            /**
             * @brief The symbolic state deviant.
             *
             */
            casadi::SX cdx;

            /**
             * @brief The actual symbolic input.
             *
             */
            casadi::SX cu;

            /**
             * @brief The generalized symbolic input (summed wrenches, used for the centroidal dynamics).
             *
             */
            casadi::SX cu_general;

            /**
             * @brief The symbolic time.
             *
             */
            casadi::SX cdt;

            /**
             * @brief The number of end effectors in the robot.
             *
             */
            int num_end_effectors_;

        private:
            /**
             * @brief The IDs that correspond to the pinocchio end effector frames.
             *
             */
            std::vector<pinocchio::FrameIndex> ee_ids_;

            /**
             * @brief Contact combination "1 2 3", where ee's 1,2, and 3 are in contact is at index with a binary value of (1110) for a system with 4 EEs.
             *
             */
            std::vector<contact::ContactCombination> contact_combinations_;

            /**
             * @brief The end effector data of the robot.
             *
             */
            contact::RobotEndEffectors ees_;

            /**
             * @brief Eigen representations of the symbolic joint variables (used for pinocchio).
             *
             */
            opt::ConfigVectorAD q_AD;

            /**
             * @brief Eigen representations of the symbolic joint velocities (used for pinocchio).
             *
             */
            opt::TangentVectorAD v_AD;
        };
    }
}