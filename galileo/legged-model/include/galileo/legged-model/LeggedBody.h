#pragma once

#include "galileo/legged-model/ContactSequence.h"
#include "galileo/math/LieAlgebra.h"
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
             * @param end_effector_names The string IDs that correspond to the pinocchio end effector frames.
             *
             */
            LeggedBody(const std::string location, const std::vector<std::string> end_effector_names);

            /**
             * @brief Provide the string IDs that correspond to the pinocchio end effector frames.
             *
             * @param ee_names The string IDs that correspond to the pinocchio end effector frames.
             */
            void setEndEffectors(const std::vector<std::string> &ee_names);

            /**
             * @brief Apply the foot velocity cost to the input cost weight matrix.
             * 
             * This cost is taken w.r.t the Jacobian at the initial configuration and acts as a hueristic to prevent large foot velocities.
             * 
             * @param R_taskspace The input cost weight matrix.
             * @param q0 The initial configuration.
             * @return Eigen::MatrixXd The updated cost weight matrix.
             */
            Eigen::MatrixXd initializeInputCostWeight(Eigen::MatrixXd R_taskspace, ConfigVector q0);

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
             * @brief Create the state derivative function to translate from the state space to the tangent space.
             *
             */
            void createFdiff();

            /**
             * @brief Create the state error function.
             *
             */
            void createErrorFunction();

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
            Model model;

            /**
             * @brief The pinocchio data of the robot.
             *
             */
            Data data;

            /**
             * @brief The symbolic model of the robot.
             *
             */
            ADModel cmodel;

            /**
             * @brief The symbolic data of the robot.
             *
             */
            ADData cdata;

            /**
             * @brief The state indices helper for the robot.
             *
             */
            std::shared_ptr<legged::LeggedRobotStates> si;

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
            ConfigVectorAD q_AD;
        };
    }
}