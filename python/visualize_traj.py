#!.venv/bin/python3

import time
import numpy as np
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

with open('python/metadata.csv', 'r') as file:
    lines = file.readlines()

# Parse the location
location = lines[0].strip().split(': ')[1][3:]  # omit "../"

# Parse the q0 array
q0_str = lines[1].strip().split(': ')[1]
q0 = np.array([float(x) for x in q0_str.split(', ')])

builder = RobotWrapper.BuildFromURDF
robot = builder(
    location,
    ["resources"],
    None,
)

robot.q0 = q0

# The pinocchio model is what we are really interested by.
model = robot.model
data = model.createData()

new_times = np.genfromtxt('python/new_times.csv', delimiter=',')
new_sol = np.genfromtxt('python/new_sol.csv', delimiter=',')
new_times = np.reshape(np.diff(new_times), (new_times.shape[0] - 1, 1))

viz = MeshcatVisualizer(robot)
viz.display(robot.q0)

def display_scene(q: np.ndarray, dt=1e-1):
    pin.framesForwardKinematics(model, data, q)
    viz.display(q)
    time.sleep(dt)

def display_traj(qs: np.ndarray, dts: np.ndarray):
    for k in range(np.size(dts, 0)):
        display_scene(qs[:, k], dts[k, 0])

while True:
    display_traj(new_sol, new_times)
    time.sleep(5)