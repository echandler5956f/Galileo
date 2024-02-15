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

            void createGeneralFunctions();

            void createGeneralDynamics();

            void createFint();

            void createFdiff();

            void fillModeDynamics();

            /**
             * @brief Get the Contact Combination object from a binary combination mask.
             *
             * @param contact_mask The binary combination mask.
             * @return contact::ContactCombination The contact combination.
             */
            contact::ContactCombination getContactCombination(int contact_mask);

            /**
             * @brief Get the End Effectors object.
             *
             */
            contact::RobotEndEffectors getEndEffectors();

            opt::Model model;

            opt::Data data;

            opt::ADModel cmodel;

            opt::ADData cdata;

            std::shared_ptr<opt::LeggedRobotStates> si;

            std::shared_ptr<contact::ContactSequence> contact_sequence;

            casadi::Function general_dynamics;

            casadi::Function fint;

            casadi::Function fdif;

            casadi::SX cx;

            casadi::SX cdx;

            casadi::SX cu;

            casadi::SX cu_general;

            casadi::SX cdt;

            /**
             * @brief The number of end effectors in the robot.
             *
             */
            int num_end_effectors_;

        private:
            casadi::SX customFint(casadi::SX x, casadi::SX dx, casadi::SX dt);

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

            opt::ConfigVectorAD q_AD;

            opt::TangentVectorAD v_AD;
        };
    }
}