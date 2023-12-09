from contacts.end_effector import EndEffector
from frozendict import frozendict
from dataclasses import dataclass


# Class for storing phase
@dataclass(frozen=True)
class Phase:
    # Dictionary which maps end effectors to contact flags (True is in contact)
    contacts: frozendict({EndEffector: bool})
    # Name of the phase (used for hashing)
    phase_name: str = None
    # If the phase has a fixed period, this value is set
    fixed_timing: float = None
    # If the phase has a variable period (determined using Casadi), this flag is True
    timing_var: bool = False
    # Number of knot points in this phase
    knot_points: int = 1
