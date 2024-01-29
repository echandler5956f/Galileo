#pragma once

#include "galileo/legged-model/EndEffector.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"

#include "galileo/opt/PhaseSequence.h"

namespace galileo
{
    namespace legged
    {
        namespace contact
        {

            /**
             * @brief A struct for holding the contact mode of the robot.
             *
             */
            struct ContactMode
            {
                /**
                 * @brief Error codes
                 *
                 */
                enum ContactModeValidity
                {
                    VALID,               // No error
                    SURFACE_NOT_DEFINED, // The surface is not defined
                    DIFFERING_SIZES      // The sizes of the contact combination and contact surfaces are different
                };
                /**
                 * @brief Maps EE's to flags (true = in contact).
                 *
                 */
                ContactCombination combination_definition;

                /**
                 * @brief Gets the contact combination for a given end effector.
                 *
                 * @param ee The end effector to get the contact combination for.
                 */
                bool &operator[](const EndEffector &ee)
                {
                    return combination_definition[ee.frame_name];
                }

                /**
                 * @brief Gets the contact combination for a given end effector.
                 *
                 * @param ee_name The name of the end effector to get the contact combination for.
                 */
                bool &operator[](const std::string &ee_name)
                {
                    return combination_definition[ee_name];
                }

                const bool &at(const EndEffector &ee) const
                {
                    return combination_definition.at(ee.frame_name);
                }

                const bool &at(const std::string &ee_name) const
                {
                    return combination_definition.at(ee_name);
                }

                /**
                 * @brief Gets the contact surface for a given end effector in this mode.
                 *
                 * @param ee_name The end effector to get the contact surface for.
                 */
                const environment::SurfaceID &getSurfaceID(const std::string &ee_name) const;

                /**
                 * @brief Gets the contact surface for a given end effector in this mode.
                 *
                 * @param ee The end effector to get the contact surface for.
                 */
                const environment::SurfaceID &getSurfaceID(const EndEffector &ee) const
                {
                    return getSurfaceID(ee.frame_name);
                }

                /**
                 * @brief Gets which surfaces the EEs are in contact with.
                 *
                 */
                std::vector<environment::SurfaceID> contact_surfaces;

                /**
                 * @brief Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE.
                 *
                 * @param validity Error code
                 */
                void MakeValid(ContactModeValidity &validity);
            };

            /**
             * @brief Class for holding simple contact sequence metadata.
             *
             */
            class ContactSequence : public PhaseSequence<ContactMode>
            {
            public:
                /**
                 * @brief Error codes.
                 *
                 */
                typedef CONTACT_SEQUENCE_ERROR PHASE_SEQUENCE_ERROR;

                /**
                 * @brief Construct a new Contact Sequence object.
                 *
                 * @param num_end_effectors The number of end effectors in the contact sequence.
                 */
                ContactSequence(int num_end_effectors) : PhaseSequence(), num_end_effectors_(num_end_effectors) {}
                
                const int &num_end_effectors() { return num_end_effectors_; }

                /**
                 * @brief Adds a phase to the contact sequence.
                 * 
                 * @param mode The contact mode of the phase.
                 * @param knot_points The number of knot points in the phase.
                 * @param dt The time step of the phase.
                */
                int addPhase(const ContactMode &mode, int knot_points, double dt) override;

                private:
                /**
                 * @brief The number of end effectors in the contact sequence.
                 */
                int num_end_effectors_;
            };
        }
    }
}