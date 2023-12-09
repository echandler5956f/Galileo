#pragma once

#include "EndEffector.h"
#include "EnvironmentSurfaces.h"

namespace acro
{
    namespace contact
    {

        struct ContactMode
        {
            enum ContactModeValidity
            {
                VALID,
                SURFACE_NOT_DEFINED,
                DIFFERING_SIZES
            };
            // Gets which EEs are in contact
            ContactCombination combination_definition;

            // Gets which surfaces the EEs are in contact with
            std::vector<environment::SurfaceID> contact_surfaces;

            // Makes the combination valid. If an EE is not in contact, it makes the corresponding contact surface NO_SURFACE
            void MakeValid(ContactModeValidity &validity);
        };

        class ContactSequence
        {
        public:
            enum CONTACT_SEQUENCE_ERROR
            {
                OK,
                NOT_IN_DT
            };

            ContactSequence(int num_end_effectors) : num_end_effectors_(num_end_effectors){}

            // Does the phase timing change? if so, then the _t0_offset and dt_ need to change.
            struct Phase
            {
                ContactMode mode;
                int knot_points = 1;
                double time_value = 1;
            };

            int addPhase(const ContactMode &mode, int knot_points, double dt);

            int getPhaseIndexAtTime(double t, CONTACT_SEQUENCE_ERROR &error_status) const;

            int getPhaseIndexAtKnot(int knot_idx, CONTACT_SEQUENCE_ERROR &error_status) const;
        
            const ContactSequence::Phase getPhase(int index) const { return phase_sequence_[index]; }

            void getPhaseAtTime(double t, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

            void getPhaseAtKnot(int knot_idx, Phase &phase, CONTACT_SEQUENCE_ERROR &error_status) const;

            int num_phases() const { return phase_sequence_.size(); }

            // we will fill this out as needed.
        // private:
            std::vector<Phase> phase_sequence_;
            std::vector<double> phase_t0_offset_;
            std::vector<int> phase_knot0_idx_;
            double dt_ = 0;
            int total_knots_ = 0;
            int num_end_effectors_;
        };
    }
}