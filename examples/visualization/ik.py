#!.venv/bin/python3

from __future__ import print_function

import numpy as np
from numpy.linalg import norm, solve
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

with open('examples/visualization/solution_data/metadata.csv', 'r') as file:
    lines = file.readlines()

# Parse the location
location = lines[0].strip().split(': ')[1].replace("../", "", 2)

builder = RobotWrapper.BuildFromURDF
robot = builder(
    location,
    ["."],
    pin.JointModelFreeFlyer(),
)

# The pinocchio model is what we are really interested by.
model = robot.model
data = model.createData()

names = model.names
for i in range(len(names)):
    print(f"{i}: {names[i]}")

q = np.zeros(model.nq)
q[2] = 0.788527 + 0.07645
print(q[2])
q[6] = 1

q[2+6] = 0

q[5 + 6] = -1.5
q[21+6] = -0.5 # hip y
q[22+6] = 1. # knee y
q[23+6] = -0.5 # ankle y

q[13+6] = 1.5
q[27+6] = -0.5 # hip y
q[28+6] = 1. # knee y
q[29+6] = -0.5 # ankle y

pin.framesForwardKinematics(model, data, q)
foot_pos = data.oMf[model.getFrameId("r_foot", pin.BODY)].translation
print(foot_pos)

viz = MeshcatVisualizer(robot)
while True:
    viz.display(q)
