#!.venv/bin/python3

from __future__ import print_function

import numpy as np
from numpy.linalg import norm, solve
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

location = "/home/quant/ros_ws/src/Galileo/resources/atlas/urdf/atlas.urdf"

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
    print(f"{i-2+7}: {names[i]}")

q = np.zeros(model.nq)
# 0, 0, 0.864977, 0, 0, 0, 1, 0, 0, 0, 0, -1.5, 0, 0, 0, 0, 0, 0, 0, 1.5, 0, 0, 0, 0, 0, 0, 0, -0.5, 1, -0.5, 0, 0, 0, -0.5, 1, -0.5, 0)
q[2] = 0.864977
q[6] = 1

q[7] = 0 # lower="-0.523599" upper="0.523599"
q[8] = 0 # lower="-0.219388" upper="0.538783"
q[9] = 0 # lower="-0.663225" upper="0.663225"

q[10] = -0.5 # lower="-1.5708" upper="0.785398"
q[11] = -1.2 # lower="-1.5708" upper="1.5708"
q[12] = 2 # lower="0" upper="3.14159"
q[13] = 0.75 # lower="0" upper="2.35619"
q[14] = 0 # lower="-3.011" upper="3.011"
q[15] = 0.25 # lower="-1.7628" upper="1.7628"
q[16] = 0 # lower="-2.9671" upper="2.9671"

q[17] = 0

q[18] = 0.5
q[19] = 1.2
q[20] = 2
q[21] = -0.75
q[22] = 0
q[23] = -0.25
q[24] = 0

q[27] = -0.5
q[28] = 1
q[29] = -0.5
q[33] = -0.5
q[34] = 1
q[35] = -0.5

# for i in range(len(q)):
#     print(q[i])

pin.framesForwardKinematics(model, data, q)
foot_pos = data.oMf[model.getFrameId("l_foot", pin.BODY)].translation
print(foot_pos)

foot_pos = data.oMf[model.getFrameId("r_foot", pin.BODY)].translation
print(foot_pos)

viz = MeshcatVisualizer(robot)
while True:
    viz.display(q)
