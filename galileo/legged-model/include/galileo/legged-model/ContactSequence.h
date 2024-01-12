#pragma once

#include "galileo/legged-model/EndEffector.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"

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
                 * @brief Gets which surfaces the EEs are in contact with.
                 *
                 */
                std::vector<environment::SurfaceID> contact_surfaces;

                /**
                 * @brief Get the surface ID for a given end effector.
                 *
                 * This function returns the surface ID for a given end effector.
                 *
                 * @param EndEffectorID The ID of the end effector.
                 * @return The surface ID for the given end effector.
                 */
                environment::SurfaceID getSurfaceIDForEE(std::string EndEffectorID) const;

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
            class ContactSequence
            {
            public:
                /**
                 * @brief Error codes.
                 *
                 */
                enum CONTACT_SEQUENCE_ERROR
                {

                    OK,       // No error
                    NOT_IN_DT // The time is out of bounds
                };

                /**
                 * @brief Construct a new Contact Sequence object.
                 *
                 * @param num_end_effectors The number of end effectors in the contact sequence.
                 */
                ContactSequence(int num_end_effectors) : num_end_effectors_(num_end_effectors) {}

                /**
                 * @brief Holds the contact mode, the number of knot points, and the time value for a certain phase. Note: Does the phase timing change? if so, then the _t0_offset and dt_ need to change.
                 *
                 */
                struct Phase
                {
                    /**
                     * @brief Contains the basic info about the current contact configuration.
                     *
                     */
                    ContactMode mode;
                    /**
                     * @brief Number of knot points for which the phase applies over.
                     *
                     */
                    int knot_points = 1;

                    /**
                     * @brief The time value for the phase.
                     */
                    double time_value = 1;
                };

                /**
                 * @brief Adds a new phase to the contact sequence.
                 *
                 * This function adds a new phase to the contact sequence with the specified contact mode,
                 * number of knot points, and time step.
                 *
                 * @param mode The contact mode of the phase.
                 * @param knot_points The number of knot points in the phase.
                 * @param dt The time step between each knot point.
                 * @return The index of the newly added phase.
                 */
                int addPhase(const ContactMode &mode, int knot_points, double dt);

                /**
                 * @brief Get the phase index at a given time.
                 *
                 * This function returns the phase index at a given time in the contact sequence.
                 *
                 * @param t The time for which to retrieve the phase index.
                 * @param error_status A reference to a CONTACT_SEQUENCE_ERROR enum that will be updated with the error status.
                 * @return The phase index at the given time.
                 */
                int getPhaseIndexAtTime(double t, CONTACT_SEQUENCE_ERROR &error_status) const;

                /**
                 * @brief Get the phase index at the specified knot index.
                 *
                 * This function returns the phase index corresponding to the given knot index.
                 *
                 * @param knot_idx The index of the knot.
                 * @param error_status The reference to the CONTACT_SEQUENCE_ERROR enum to store the error status.
                 * @return The phase index at the specified knot index.
                 */
                int getPhaseIndexAtKnot(int knot_idx, CONTACT_SEQUENCE_ERROR &error_status) const;

                /**
                 * @brief Get the phase at the specified index.
                 *
                 * @param index The index of the phase to retrieve.
                 * @return The phase at the specified index.
                 */
                const ContactSequence::Phase getPhase(int index) const { return phase_sequence_[index]; }

                /**
                 * @brief Get the phase at a given time.
                 *
                 * This function retrieves the phase at a specified time.
                 *
                 * @param t The time at which to retrieve the phase.
                 * @param phase [out] The phase at the specified time.
                 * @param error_status [out] The error status of the operation.
                 */
                void getPhaseAtTime(double t, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

                /**
                 * @brief Get the phase at a specific knot index.
                 *
                 * This function retrieves the phase at the specified knot index.
                 *
                 * @param knot_idx The index of the knot.
                 * @param phase [out] The phase at the specified knot index.
                 * @param error_status [out] The error status of the operation.
                 */
                void getPhaseAtKnot(int knot_idx, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

                /**
                 * @brief Get the number of phases in the contact sequence.
                 *
                 * @return The number of phases.
                 */
                int num_phases() const { return phase_sequence_.size(); }

                /**
                 * @brief A vector of Phase objects.
                 *
                 * This vector represents a sequence of phases.
                 */
                std::vector<Phase> phase_sequence_;

                /**
                 * @brief Struct representing the global phase offset.
                 *
                 * This struct contains the time offset (t0_offset) and knot offset (knot0_offset)
                 * for a contact sequence.
                 */
                struct GlobalPhaseOffset
                {
                    double t0_offset; /**< Time offset */
                    int knot0_offset; /**< Knot offset */
                };

                /**
                 * @brief A vector of GlobalPhaseOffset objects.
                 *
                 * This vector stores GlobalPhaseOffset objects, which represent phase offsets in a contact sequence.
                 */
                std::vector<GlobalPhaseOffset> phase_offset_;

                /**
                 * @brief The time step used for the contact sequence.
                 */
                double dt_ = 0;

                /**
                 * @brief The total number of knots.
                 */
                int total_knots_ = 0;

                /**
                 * @brief The number of end effectors in the contact sequence.
                 */
                int num_end_effectors_;
            };
        }
    }
}