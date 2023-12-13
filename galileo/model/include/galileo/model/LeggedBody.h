#pragma once

#include "galileo/model/EndEffector.h"
#include <pinocchio/multibody/model.hpp>
#include <vector>
#include <bitset>

namespace galileo
{
    namespace model
    {
        class LeggedBody : public pinocchio::Model
        {
        public:
            /**
             * @brief Construct a new Legged Body object
             *
             */
            LeggedBody() : pinocchio::Model() {}

            /**
             * @brief Provide the string IDs that correspond to the pinocchio end effector frames.
             *
             * @param ee_names
             */
            void setEndEffectors(const std::vector<std::string> &ee_names);

            /**
             * @brief Generate combinations of contacts.
             *
             */
            void GenerateContactCombination();

            /**
             * @brief Get the Contact Combination object from a binary combination mask
             *
             * @param contact_mask
             * @return contact::ContactCombination
             */
            contact::ContactCombination getContactCombination(int contact_mask) { return contact_combinations_[contact_mask]; }

        private:
            /**
             * @brief
             *
             */
            std::vector<std::string> ee_names_;

            /**
             * @brief Contact combination "1 2 3", where ee's 1,2, and 3 are in contact is at index with a binary value of (1110) for a system with 4 EEs
             *
             */
            std::vector<contact::ContactCombination> contact_combinations_;

            /**
             * @brief
             *
             */
            contact::RobotEndEffectors ees_;

            /**
             * @brief
             *
             */
            int num_end_effectors_;
        };
    }
}