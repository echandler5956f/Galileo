#!.venv/bin/python3

import numpy as np
import casadi as ca
from frozendict import frozendict
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

from contacts import *
from trajopt import *


builder = RobotWrapper.BuildFromURDF
robot = builder(
    "resources/urdf/huron_cheat.urdf",
    ["resources"],
    None,
)

robot.q0 = np.array(
    [
        0,
        0,
        1.0627,
        0,
        0,
        0,
        1,
        0.0000,
        0.0000,
        -0.3207,
        0.7572,
        -0.4365,
        0.0000,
        0.0000,
        0.0000,
        -0.3207,
        0.7572,
        -0.4365,
        0.0000,
    ]
)

# The pinocchio model is what we are really interested by.
model = robot.model
data = model.createData()

viz = MeshcatVisualizer(robot)
viz.display(robot.q0)

params = Parameters(robot.q0)

l_name = "l_foot_v_ft_link"
r_name = "r_foot_v_ft_link"

# Define two end effectors
ee_left_foot = EndEffector(
    frame_name=l_name, frame_id=model.getFrameId(l_name), type_6D=True
)
ee_right_foot = EndEffector(
    frame_name=r_name, frame_id=model.getFrameId(r_name), type_6D=True
)
phase_1_period = 0.25
phase_1_knot_points = 25
phase_1_fixed_timing = phase_1_period / phase_1_knot_points
contacts_phase1 = frozendict({ee_left_foot: True, ee_right_foot: True})
phase1 = Phase(
    contacts=contacts_phase1,
    phase_name="phase_1",
    fixed_timing=phase_1_fixed_timing,
    knot_points=phase_1_knot_points,
)
phase_2_period = 0.25
phase_2_knot_points = 25
phase_2_fixed_timing = phase_2_period / phase_2_knot_points
contacts_phase2 = frozendict({ee_left_foot: True, ee_right_foot: False})
phase2 = Phase(
    contacts=contacts_phase2,
    phase_name="phase_2",
    fixed_timing=phase_2_fixed_timing,
    knot_points=phase_2_knot_points,
)
phase_3_period = 0.05
phase_3_knot_points = 5
phase_3_fixed_timing = phase_3_period / phase_3_knot_points
contacts_phase3 = frozendict({ee_left_foot: True, ee_right_foot: True})
phase3 = Phase(
    contacts=contacts_phase3,
    phase_name="phase_3",
    fixed_timing=phase_3_fixed_timing,
    knot_points=phase_3_knot_points,
)

# Create a contact sequence and add phases
contact_seq = ContactSequence()
contact_seq.add_phase(phase1)
contact_seq.add_phase(phase2)
# contact_seq.add_phase(phase3)

traj_opt = CentroidalTrajOpt(model, data, viz, params, contact_seq)
traj_opt.compute_casadi_graphs()
print("Finished computing casadi graphs")
traj_opt.solve_problem()
