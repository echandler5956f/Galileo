#!.venv/bin/python3

import time
import numpy as np
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

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

new_times = np.genfromtxt('python/new_times.csv', delimiter=',')
new_sol = np.genfromtxt('python/new_sol.csv', delimiter=',')
new_times = np.reshape(np.diff(new_times), (new_times.shape[0] - 1, 1))

viz = MeshcatVisualizer(robot)
viz.display(robot.q0)

ee = ["r_foot_v_ft_link", "l_foot_v_ft_link"]

def display_scene(q: np.ndarray, dt=1e-1):
    """
    Given the robot configuration, display:
    - the robot
    - a box representing endEffector_ID
    - a box representing Mtarget
    """
    pin.framesForwardKinematics(model, data, q)
    viz.applyConfiguration(ee[0], data.oMf[model.getFrameId(ee[0])])
    viz.applyConfiguration(ee[1], data.oMf[model.getFrameId(ee[1])])
    viz.display(q)
    time.sleep(dt)

def display_traj(qs: np.ndarray, dts: np.ndarray):
    for k in range(np.size(dts, 0)):
        display_scene(qs[:, k], dts[k, 0])

while True:
    display_traj(new_sol, new_times)
    time.sleep(5)