#pragma once

#include "galileo/legged-model/EndEffector.h"
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/algorithm/joint-configuration.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <Eigen/LU>
#include <pinocchio/parsers/urdf.hpp>
#include <vector>
#include <bitset>

namespace galileo
{
    namespace legged
    {
        /**
         * @brief A class for holding the robot model.
         *
         */
        template <class Scalar>
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
            LeggedBody(const std::string location, const int num_ees, const std::string end_effector_names[])
            {
                this->model = pinocchio::ModelTpl<Scalar>();
                std::vector<std::string> ee_name_vect;
                ee_name_vect.resize(num_ees);
                for (int i = 0; i < num_ees; ++i)
                {
                    ee_name_vect[i] = end_effector_names[i];
                }
                pinocchio::urdf::buildModel(location, this->model);
                this->setEndEffectors(ee_name_vect);
                this->GenerateContactCombination();
            }

            /**
             * @brief Provide the string IDs that correspond to the pinocchio end effector frames.
             *
             * @param ee_names The string IDs that correspond to the pinocchio end effector frames.
             */
            void setEndEffectors(const std::vector<std::string> &ee_names)
            {
                auto data = pinocchio::DataTpl<Scalar>(this->model);
                for (int i = 0; i < ee_names.size(); ++i)
                {
                    auto &ee_name = ee_names[i];
                    assert(this->model.existFrame(ee_name));
                    // Throw an error otherwise.
                    std::shared_ptr<contact::EndEffector> ee_obj_ptr = std::make_shared<contact::EndEffector>();
                    ee_obj_ptr->frame_name = ee_name;
                    ee_obj_ptr->frame_id = this->model.getFrameId(ee_name);

                    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix(6, this->model.nv);
                    matrix.setZero();
                    pinocchio::computeFrameJacobian(model, data, pinocchio::neutral(model), ee_obj_ptr->frame_id, pinocchio::LOCAL, matrix);
                    Eigen::FullPivLU<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> lu_decomp(matrix);
                    ee_obj_ptr->is_6d = lu_decomp.rank() == 6;
                    ee_obj_ptr->local_ee_idx = i;

                    ees_.insert({ee_name, ee_obj_ptr});
                }
                num_end_effectors_ = ee_names.size();
                ee_names_ = ee_names;
            }

            /**
             * @brief Generate combinations of contacts.
             *
             */
            void GenerateContactCombination()
            {
                std::vector<contact::ContactCombination> contact_combinations;
                // Generate the "basic" (no contact) contact combination.
                contact::ContactCombination basic_cc;
                for (auto &ee_name : ee_names_)
                    basic_cc.insert({ee_name, false});

                // The power set can be gotten by all the binary values between 0 and 2^n-1.
                // That is, 0000, 0001, 0010, 0011, ... for 4 EEs.
                int num_combinations = pow(2, num_end_effectors_);

                contact_combinations.resize(num_combinations);
                for (uint binary_value_combination = 0; binary_value_combination < num_combinations; binary_value_combination++)
                {
                    // Copy the no contact cc
                    contact::ContactCombination new_contact_combination = basic_cc;
                    // And set the value for each ee in contact to true
                    std::bitset<32> bit_mask = 1;
                    std::bitset<32> binary_value_combination_bits = binary_value_combination;

                    for (int i = 0; i < num_end_effectors_; i++)
                    {
                        bool ee_i_is_in_contact = (bit_mask & binary_value_combination_bits).any();
                        // bit shift
                        bit_mask <<= 1;
                        new_contact_combination[ee_names_[i]] = ee_i_is_in_contact;
                    }
                    contact_combinations[binary_value_combination] = new_contact_combination;
                }
                contact_combinations_ = contact_combinations;
            }

            /**
             * @brief Get the Contact Combination object from a binary combination mask.
             *
             * @param contact_mask The binary combination mask.
             * @return contact::ContactCombination The contact combination.
             */
            contact::ContactCombination getContactCombination(int contact_mask) { return contact_combinations_[contact_mask]; }

            /**
             * @brief Get the End Effectors object.
             *
             */
            contact::RobotEndEffectors getEndEffectors() { return ees_; }

            /**
             * @brief The robot model.
             *
             */
            pinocchio::ModelTpl<Scalar> model;

        private:
            /**
             * @brief The string IDs that correspond to the pinocchio end effector frames.
             *
             */
            std::vector<std::string> ee_names_;

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
             * @brief The number of end effectors in the robot.
             *
             */
            int num_end_effectors_;
        };
    }
}