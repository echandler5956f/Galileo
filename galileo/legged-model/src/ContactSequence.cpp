#include "galileo/legged-model/ContactSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace contact
        {
            const environment::SurfaceID &ContactMode::getSurfaceID(const std::string &ee_name) const
            {
                const auto &map_element = combination_definition.find(ee_name);

                if (map_element == combination_definition.end())
                {
                    throw std::runtime_error(std::string("'End Effector " + ee_name + " not defined in contact mode!'"));
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

        int ContactSequence::addPhase(const ContactMode &mode, int knot_points, double dt)
        {
            auto validity_mode = mode;
            ContactMode::ContactModeValidity validity;

            validity_mode.MakeValid(validity);

            if (validity != ContactMode::VALID)
            {
                throw std::runtime_error(std::string("'Contact is not valid!'"));
            }
            
            PhaseSequence<ContactMode>::addPhase(validity_mode, knot_points, dt);
        }

        }
    }
}