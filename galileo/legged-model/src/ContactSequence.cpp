#include "galileo/legged-model/ContactSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace contact
        {
            const environment::SurfaceID &ContactMode::getSurfaceID(pinocchio::FrameIndex ee_id) const
            {
                const auto &map_element = combination_definition.find(ee_id);

                if (map_element == combination_definition.end())
                {
                    throw std::runtime_error(std::string("End Effector " + std::to_string(ee_id) + " not defined in contact mode!"));
                }
                uint index = std::distance(combination_definition.begin(), map_element);

                return contact_surfaces[index];
            }

            // Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE
            void ContactMode::MakeValid(ContactModeValidity &validity)
            {
                if (combination_definition.size() != contact_surfaces.size())
                {
                    validity = ContactModeValidity::DIFFERING_SIZES;
                    return;
                }

                ContactCombination::iterator it = combination_definition.begin();

                for (uint i = 0; i < combination_definition.size(); i++)
                {
                    bool is_in_contact = (*it).second;

                    if (!is_in_contact)
                    {
                        // If it is not in contact, make the contact surface NO_SURFACE
                        contact_surfaces[i] = environment::NO_SURFACE;
                    }

                    bool surface_is_not_defined = (contact_surfaces[i] == environment::NO_SURFACE);

                    if (is_in_contact && surface_is_not_defined)
                    {
                        validity = ContactModeValidity::SURFACE_NOT_DEFINED;
                        return;
                    }

                    if (i < combination_definition.size() - 1)
                        it++;
                }

                validity = ContactMode::ContactModeValidity::VALID;
            }

            int ContactMode::numEndEffectorsInContact() const
            {
                int count = 0;
                for (const auto &pair : combination_definition)
                {
                    if (pair.second) // if the end effector is in contact
                    {
                        count++;
                    }
                }
                return count;
            }

            int ContactSequence::addPhase(const ContactMode &mode, int knot_points, double dt)
            {
                auto validity_mode = mode;
                ContactMode::ContactModeValidity validity;

                validity_mode.MakeValid(validity);

                if (validity != ContactMode::VALID)
                {
                    throw std::runtime_error(std::string("Contact is not valid!"));
                }
                return commonAddPhase(mode, knot_points, dt);
            }

            int ContactSequence::numEndEffectorsInContactAtPhase(int phase_index) const
            {
                if (phase_index < 0 || size_t(phase_index) >= phase_sequence_.size())
                {
                    throw std::out_of_range("Phase index out of range");
                }

                return phase_sequence_[phase_index].mode.numEndEffectorsInContact();
            }
        }
    }
}