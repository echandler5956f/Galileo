from contacts.end_effector import EndEffector
from contacts.phase import Phase
from collections import defaultdict
from typing import List


# Class for operating on contacts
class ContactSequence:
    def __init__(self):
        self.sequence: List[Phase] = []
        self.cumulative_knots: List[int] = []

    # Add a new phase to the sequence
    # Cumulative_knots[-1] returns the total number of knot points in the sequence
    def add_phase(self, phase: Phase):
        self.sequence.append(phase)

        previous_knots = self.cumulative_knots[-1] if self.cumulative_knots else 0
        self.cumulative_knots.append(previous_knots + phase.knot_points)

    # Get a list of the end effectors (useful for iterating over)
    def get_all_end_effectors(self) -> List[EndEffector]:
        end_effectors = set()
        for phase in self.sequence:
            end_effectors.update(phase.contacts.keys())
        return list(end_effectors)

    # Get the number of end effectors in contact during a phase
    def get_num_stance_during_phase(self, phase_idx: int) -> int:
        count = 0
        for flag in self.sequence[phase_idx].contacts.values():
            count += int(flag)
        return count

    # Get phase from a knot point index
    def get_phase(self, k: int) -> Phase:
        return self.sequence[self.get_phase_idx(k)]

    # Return whether an end effector is in contact at a certain knot point
    def is_in_contact(self, end_effector: EndEffector, k: int) -> bool:
        phase = self.get_phase(k)
        return phase.contacts.get(end_effector, False)

    def is_new_contact(self, k: int, end_effector: EndEffector) -> bool:
        i = self.get_phase_idx(k)
        if self.is_in_contact(end_effector, k) and i > 0:
            if (
                not self.sequence[i - 1].contacts[end_effector]
                and k - self.cumulative_knots[i - 1] == 0
            ):
                return True
        elif i == 0:
            return True
        return False

    # Get the frame id from an end effector object (used with pinocchio)
    def get_frame_id(self, end_effector: EndEffector) -> int:
        return end_effector.frame_id

    # Get the flag which specifies whether the end effector has 6 DOF (True) or 3 DOF
    def get_contact_type(self, end_effector: EndEffector) -> bool:
        return end_effector.type_6D

    # Convert the contact type flag into a number (True = 6, False = 3)
    # Used when determining input variable size with casadi
    def get_contact_size(self, end_effector: EndEffector) -> int:
        return 3 * int(end_effector.type_6D) + 3

    # Get the total number of phases in the sequence
    def get_num_phases(self) -> int:
        return len(self.sequence)

    # Convert a knot point into a phase index
    def get_phase_idx(self, k: int) -> int:
        for idx, cum_knot in enumerate(self.cumulative_knots):
            if k < cum_knot:
                return idx
        return None

    def get_time_vec(self, dts: defaultdict(list)):
        t = defaultdict(list)
        k = 0
        tc = 0.0
        for i in range(self.get_num_phases()):
            phase = self.sequence[i]
            t[phase.phase_name] = []
            for ki in range(phase.knot_points):
                t[phase.phase_name].append(tc + dts[i])
                tc += dts[i]
                k += 1
        return t
