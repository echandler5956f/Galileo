from dataclasses import dataclass


# Class for storing end efector information
@dataclass(frozen=True)
class EndEffector:
    # Name of the frame in the URDF
    frame_name: str
    # ID of the frame in pinocchio
    frame_id: int
    # True if the end effector is 6-DOF, False if 3-DOF
    type_6D: bool
