#include "LeggedBody.h"

// Provide the string IDs that correspond to the pinocchio end effector frames.
void acro::model::LeggedBody::setEndEffectors(const std::vector<std::string> &ee_names)
{
    for (int i = 0; i < ee_names.size(); ++i)
    {
        auto &ee_name = ee_names[i];
        assert(existFrame(ee_name));
        // Throw an error otherwise.
        std::shared_ptr<contact::EndEffector> ee_obj_ptr = std::make_shared<contact::EndEffector>();
        ee_obj_ptr->frame_name = ee_name;
        ee_obj_ptr->frame_id = getFrameId(ee_name);
        
        // todo : use the jacobian from the body frame to find the actual DOFs of this frame.
        ee_obj_ptr->is_6d = false;

        ee_obj_ptr->local_ee_idx = i;

        ees_.insert({ee_name, ee_obj_ptr});
    }
    num_end_effectors_ = ee_names.size();
    ee_names_ = ee_names;
}

// Generate combinations of contacts.
void acro::model::LeggedBody::GenerateContactCombination()
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
            //bit shift
            bit_mask<<=1;
            new_contact_combination[ee_names_[i]] = ee_i_is_in_contact;
        }
        contact_combinations[binary_value_combination] = new_contact_combination;
    }
    contact_combinations_ = contact_combinations;
}