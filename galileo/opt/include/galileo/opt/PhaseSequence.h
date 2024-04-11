#pragma once

#include <vector>
#include <pinocchio/autodiff/casadi.hpp>

namespace galileo
{
    namespace opt
    {
        /**
         * @brief Class for holding simple Phase sequence metadata.
         *
         */
        template <typename MODE_T>
        class PhaseSequence
        {
        public:
            /**
             * @brief Error codes.
             *
             */
            enum PHASE_SEQUENCE_ERROR
            {
                OK,       // No error
                NOT_IN_DT // The time is out of bounds
            };

            struct Phase;

            /**
             * @brief Construct a new Phase Sequence object.
             *
             */
            PhaseSequence() {}

            PhaseSequence(const PhaseSequence &other);

            PhaseSequence(const std::vector<MODE_T> &modes, const std::vector<int> &knot_points, const std::vector<double> &dt);

            PhaseSequence(const std::vector<Phase> &phases);

            PhaseSequence &operator=(const PhaseSequence &other);

            PhaseSequence &operator=(const std::vector<Phase> &phases);

            /**
             * @brief Destroy the Phase Sequence object.
             *
             */
            ~PhaseSequence() = default;

            /**
             * @brief Holds the mode, the number of knot points, and the time value for a certain phase. Note: Does the phase timing change? if so, then the _t0_offset and dt_ need to change.
             *
             */
            struct Phase
            {
                /**
                 * @brief Contains the basic info about the current mode.
                 *
                 */
                MODE_T mode;

                /**
                 * @brief Dynamics of this phase.
                 *
                 */
                casadi::Function phase_dynamics;

                /**
                 * @brief Cost function for this phase.
                 *
                 */
                casadi::Function phase_cost;

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
             * @brief Get the phase index at a given time.
             *
             * This function returns the phase index at a given time in the sequence.
             *
             * @param t The time for which to retrieve the phase index.
             * @param error_status A reference to a PHASE_SEQUENCE_ERROR enum that will be updated with the error status.
             * @return The phase index at the specified time.
             */
            int getPhaseIndexAtTime(double t, PHASE_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the phase index at the specified knot index.
             *
             * This function returns the phase index corresponding to the given knot index.
             *
             * @param knot_idx The index of the knot.
             * @param error_status The reference to the PHASE_SEQUENCE_ERROR enum to store the error status.
             * @return The phase index at the specified knot index.
             */
            int getPhaseIndexAtKnot(int knot_idx, PHASE_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the phase at the specified index.
             *
             * @param index The index of the phase to retrieve.
             * @return The phase at the specified index.
             */
            const Phase getPhase(int index) const { return phase_sequence_[index]; }

            /**
             * @brief Get the phase at a given time.
             *
             * This function retrieves the phase at a specified time.
             *
             * @param t The time at which to retrieve the phase.
             * @param phase [out] The phase at the specified time.
             * @param error_status [out] The error status of the operation.
             */
            void getPhaseAtTime(double t, Phase &phase, PHASE_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the time at a given phase index.
             *
             * This function retrieves the time at a specified phase index.
             *
             * @param phase_idx The index of the phase.
             * @param t [out] The time at the specified phase index.
             * @param error_status [out] The error status of the operation.
             */
            void getTimeAtPhase(int phase_idx, double &t, PHASE_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the phase at a specific knot index.
             *
             * This function retrieves the phase at the specified knot index.
             *
             * @param knot_idx The index of the knot.
             * @param phase [out] The phase at the specified knot index.
             * @param error_status [out] The error status of the operation.
             */
            void getPhaseAtKnot(int knot_idx, Phase &phase, PHASE_SEQUENCE_ERROR &error_status) const;

            /**
             * @brief Get the number of phases in the Phase sequence.
             *
             * @return The number of phases.
             */
            int getNumPhases() const { return phase_sequence_.size(); }

            /**
             * @brief Get the time step used for the Phase sequence.
             *
             * @return const double& The time step.
             */
            const double &getDT() { return dt_; }

            /**
             * @brief Get the total number of knots in the Phase sequence.
             *
             * @return const int& The total number of knots.
             */
            const int &getTotalKnots() { return total_knots_; }

            /**
             * @brief Get the Phase sequence.
             *
             * @return const std::vector<Phase>& The Phase sequence.
             */
            const std::vector<Phase> &getPhaseSequence() { return phase_sequence_; }

            /**
             * @brief Cumulative sum of the phase timings.
             *
             * @return const std::vector<double> The cumulative sum of the phase timings.
             */
            const std::vector<double> getPhaseStartTimes() const;

            /**
             * @brief Get the phases in the Phase sequence.
             *
             */
            const std::vector<Phase> &getPhases() const { return phase_sequence_; }

            /**
             * @brief Adds dynamics to the specified phase.
             */
            const Phase &FillPhaseDynamics(int phase_idx, const casadi::Function &phase_dynamics)
            {
                phase_sequence_[phase_idx].phase_dynamics = phase_dynamics;
                return phase_sequence_[phase_idx];
            }

            /**
             * @brief Adds a phase cost to the specified phase.
             */
            const Phase &FillPhaseCost(int phase_idx, const casadi::Function &phase_cost)
            {
                phase_sequence_[phase_idx].phase_cost = phase_cost;
                return phase_sequence_[phase_idx];
            }

        protected:
            /**
             * @brief A vector of Phase objects.
             *
             * This vector represents a sequence of phases.
             */
            std::vector<Phase> phase_sequence_;

            /**
             * @brief Adds a new phase to the Phase sequence.
             *
             * This function adds a new phase to the sequence with the specified mode,
             * number of knot points, and time step.
             *
             * @param mode The mode of the phase.
             * @param knot_points The number of knot points in the phase.
             * @param dt The time step between each knot point.
             * @return The index of the newly added phase.
             */
            int commonAddPhase(const MODE_T &mode, int knot_points, double dt);

            /**
             * @brief Adds a new phase to the Phase sequence.
             *
             * This function adds a new phase to the sequence with the specified mode,
             * number of knot points, and time step.
             *
             * @param mode The mode of the phase.
             * @param knot_points The number of knot points in the phase.
             * @param dt The time step between each knot point.
             * @param phase_dynamics The dynamics of the phase. If you do not set the dynamics here, you must do so later manually!
             * @return The index of the newly added phase.
             */
            int commonAddPhase(const MODE_T &mode, int knot_points, double dt, casadi::Function phase_dynamics);

            /**
             * @brief Struct representing the global phase offset.
             *
             * This struct contains the time offset (t0_offset) and knot offset (knot0_offset)
             * for a Phase sequence.
             */
            struct GlobalPhaseOffset
            {
                /**
                 * @brief Time offset.
                 *
                 */
                double t0_offset;

                /**
                 * @brief Knot offset.
                 *
                 */
                int knot0_offset;
            };

            /**
             * @brief A vector of GlobalPhaseOffset objects.
             *
             * This vector stores GlobalPhaseOffset objects, which represent phase offsets in a Phase sequence.
             */
            std::vector<GlobalPhaseOffset> phase_offset_;

            /**
             * @brief The time step used for the Phase sequence.
             */
            double dt_ = 0;

            /**
             * @brief The total number of knots.
             */
            int total_knots_ = 0;
        };

        template <typename MODE_T>
        PhaseSequence<MODE_T>::PhaseSequence(const PhaseSequence &other)
        {
            phase_sequence_ = other.phase_sequence_;
            phase_offset_ = other.phase_offset_;
            dt_ = other.dt_;
            total_knots_ = other.total_knots_;
        }

        template <typename MODE_T>
        PhaseSequence<MODE_T>::PhaseSequence(const std::vector<MODE_T> &modes, const std::vector<int> &knot_points, const std::vector<double> &dt)
        {
            for (int i = 0; i < modes.size(); i++)
            {
                commonAddPhase(modes[i], knot_points[i], dt[i]);
            }
        }

        template <typename MODE_T>
        PhaseSequence<MODE_T>::PhaseSequence(const std::vector<Phase> &phases)
        {
            for (int i = 0; i < phases.size(); i++)
            {
                commonAddPhase(phases[i].mode, phases[i].knot_points, phases[i].time_value);
            }
        }

        template <typename MODE_T>
        PhaseSequence<MODE_T> &PhaseSequence<MODE_T>::operator=(const PhaseSequence &other)
        {
            phase_sequence_ = other.phase_sequence_;
            phase_offset_ = other.phase_offset_;
            dt_ = other.dt_;
            total_knots_ = other.total_knots_;
            return *this;
        }

        template <typename MODE_T>
        PhaseSequence<MODE_T> &PhaseSequence<MODE_T>::operator=(const std::vector<Phase> &phases)
        {
            for (int i = 0; i < phases.size(); i++)
            {
                commonAddPhase(phases[i].mode, phases[i].knot_points, phases[i].time_value);
            }
            return *this;
        }

        template <typename MODE_T>
        int PhaseSequence<MODE_T>::commonAddPhase(const MODE_T &mode, int knot_points, double dt)
        {
            Phase new_phase;
            GlobalPhaseOffset new_phase_offset;
            new_phase.mode = mode;

            new_phase.knot_points = knot_points;
            new_phase.time_value = dt;

            new_phase_offset.t0_offset = dt_;
            new_phase_offset.knot0_offset = total_knots_;

            phase_sequence_.push_back(new_phase);
            phase_offset_.push_back(new_phase_offset);
            dt_ += dt;
            total_knots_ += knot_points;
            return phase_sequence_.size() - 1;
        }

        template <typename MODE_T>
        int PhaseSequence<MODE_T>::commonAddPhase(const MODE_T &mode, int knot_points, double dt, casadi::Function phase_dynamics)
        {
            Phase new_phase;
            GlobalPhaseOffset new_phase_offset;
            new_phase.mode = mode;

            new_phase.knot_points = knot_points;
            new_phase.time_value = dt;
            new_phase.phase_dynamics = phase_dynamics;

            new_phase_offset.t0_offset = dt_;
            new_phase_offset.knot0_offset = total_knots_;

            phase_sequence_.push_back(new_phase);
            phase_offset_.push_back(new_phase_offset);
            dt_ += dt;
            total_knots_ += knot_points;
            return phase_sequence_.size() - 1;
        }

        // TODO: There's a very strange bug with this function! Accessing dt_ causes a segfault when called in LeggedRos, even though it is properly initialized....
        template <typename MODE_T>
        int PhaseSequence<MODE_T>::getPhaseIndexAtTime(double t, PHASE_SEQUENCE_ERROR &error_status) const
        {
            if ((t < 0) || (t > dt_))
            {
                error_status = PHASE_SEQUENCE_ERROR::NOT_IN_DT;
                return -1;
            }

            for (int i = getNumPhases() - 1; i > 0; i--)
            {
                bool is_in_phase_i = (t >= phase_offset_[i].t0_offset);
                std::cout << "t: " << t << " t0_offset: " << phase_offset_[i].t0_offset << " is_in_phase_i: " << is_in_phase_i << std::endl;
                if (is_in_phase_i)
                {
                    error_status = PHASE_SEQUENCE_ERROR::OK;
                    return i;
                }
            }
            return -1;
        }

        template <typename MODE_T>
        int PhaseSequence<MODE_T>::getPhaseIndexAtKnot(int knot_idx, PHASE_SEQUENCE_ERROR &error_status) const
        {
            if ((knot_idx < 0) || (knot_idx >= total_knots_))
            {
                error_status = PHASE_SEQUENCE_ERROR::NOT_IN_DT;
                return -1;
            }

            for (int i = getNumPhases() - 1; i > 0; i--)
            {
                bool is_in_phase_i = (knot_idx >= phase_offset_[i].knot0_offset);
                if (is_in_phase_i)
                {
                    error_status = PHASE_SEQUENCE_ERROR::OK;
                    return i;
                }
            }
            return -1;
        }

        // TODO: There's a very strang bug in "getPhaseIndexAtTime"... See above.
        template <typename MODE_T>
        void PhaseSequence<MODE_T>::getPhaseAtTime(double t, Phase &phase, PHASE_SEQUENCE_ERROR &error_status) const
        {
            int phase_index = getPhaseIndexAtTime(t, error_status);
            if (error_status != PHASE_SEQUENCE_ERROR::OK)
            {
                return;
            }
            phase = phase_sequence_[phase_index];
        }

        template <typename MODE_T>
        void PhaseSequence<MODE_T>::getTimeAtPhase(int phase_idx, double &t, PHASE_SEQUENCE_ERROR &error_status) const
        {
            if ((phase_idx < 0) || (phase_idx >= getNumPhases()))
            {
                error_status = PHASE_SEQUENCE_ERROR::NOT_IN_DT;
                return;
            }

            t = phase_offset_[phase_idx].t0_offset;
            error_status = PHASE_SEQUENCE_ERROR::OK;
        }

        template <typename MODE_T>
        void PhaseSequence<MODE_T>::getPhaseAtKnot(int knot_idx, Phase &phase, PHASE_SEQUENCE_ERROR &error_status) const
        {
            int phase_index = getPhaseIndexAtKnot(knot_idx, error_status);
            if (error_status != PHASE_SEQUENCE_ERROR::OK)
            {
                return;
            }

            phase = phase_sequence_[phase_index];
        }

        template <typename MODE_T>
        const std::vector<double> PhaseSequence<MODE_T>::getPhaseStartTimes() const
        {
            std::vector<double> phase_timing;
            phase_timing.push_back(0);
            PHASE_SEQUENCE_ERROR error_status;
            for (int i = 0; i < getNumPhases(); i++)
            {
                double t_i;
                getTimeAtPhase(i, t_i, error_status);
                phase_timing.push_back(t_i);
            }
            return phase_timing;
        }
    }
}