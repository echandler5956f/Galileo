#pragma once

#include "EndEffector.h"
#include "EnvironmentSurfaces.h"

namespace acro
{
    namespace contact
    {

        /**
         * @brief
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
                VALID,
                SURFACE_NOT_DEFINED,
                DIFFERING_SIZES
            };
            /**
             * @brief Maps EE's to flags (true = in contact)
             *
             */
            ContactCombination combination_definition;

            /**
             * @brief Gets which surfaces the EEs are in contact with
             *
             */
            std::vector<environment::SurfaceID> contact_surfaces;

            /**
             * @brief Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE
             *
             * @param validity Error code
             */
            void MakeValid(ContactModeValidity &validity);
        };

        /**
         * @brief Class for holding simple contact sequence metadata
         *
         */
        class ContactSequence
        {
        public:
            enum CONTACT_SEQUENCE_ERROR
            {
                OK,
                NOT_IN_DT
            };

            ContactSequence(int num_end_effectors) : num_end_effectors_(num_end_effectors) {}

            /**
             * @brief Holds the contact mode, the number of knot points, and the time value for a certain phase. Note: Does the phase timing change? if so, then the _t0_offset and dt_ need to change.
             *
             */
            struct Phase
            {
                /**
                 * @brief Contains the basic info about the current contact configuration
                 *
                 */
                ContactMode mode;
                /**
                 * @brief Number of knot points for which the phase applies over
                 *
                 */
                int knot_points = 1;
                /**
                 * @brief Duration of the phase
                 *
                 */
                double time_value = 1;
            };

            /**
             * @brief
             *
             * @param mode
             * @param knot_points
             * @param dt
             * @return int
             */
            int addPhase(const ContactMode &mode, int knot_points, double dt);

            /**
             * @brief Get the Phase Index At Time object
             *
             * @param t
             * @param error_status
             * @return int
             */
            int getPhaseIndexAtTime(double t, CONTACT_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the Phase Index At Knot object
             *
             * @param knot_idx
             * @param error_status
             * @return int
             */
            int getPhaseIndexAtKnot(int knot_idx, CONTACT_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the Phase object
             *
             * @param index
             * @return const ContactSequence::Phase
             */
            const ContactSequence::Phase getPhase(int index) const { return phase_sequence_[index]; }

            /**
             * @brief Get the Phase At Time object
             *
             * @param t
             * @param phase
             * @param error_status
             */
            void getPhaseAtTime(double t, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the Phase At Knot object
             *
             * @param knot_idx
             * @param phase
             * @param error_status
             */
            void getPhaseAtKnot(int knot_idx, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief
             *
             * @return int
             */
            int num_phases() const { return phase_sequence_.size(); }

            // we will fill this out as needed.
            // private:
            /**
             * @brief
             *
             */
            std::vector<Phase> phase_sequence_;

            /**
             * @brief
             *
             */
            std::vector<double> phase_t0_offset_;

            /**
             * @brief
             *
             */
            std::vector<int> phase_knot0_idx_;

            /**
             * @brief
             *
             */
            double dt_ = 0;

            /**
             * @brief
             *
             */
            int total_knots_ = 0;

            /**
             * @brief
             *
             */
            int num_end_effectors_;
        };
    }
}