#!.venv/bin/python3

from __future__ import print_function

import numpy as np
from numpy.linalg import norm, solve
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

location = "/home/quant/ros_ws/src/Galileo/resources/huron2/urdf/huron2.urdf"

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
q[2] = 1.1298345
q[6] = 1

q[9] = 0.3207
q[10] = -0.7572
q[11] = 0.4365
q[15] = 0.3207
q[16] = -0.7572
q[17] = -0.4365

print(q)

pin.framesForwardKinematics(model, data, q)
foot_pos = data.oMf[model.getFrameId("l_foot_v_ft_link", pin.BODY)].translation
print(foot_pos)

foot_pos = data.oMf[model.getFrameId("r_foot_v_ft_link", pin.BODY)].translation
print(foot_pos)

viz = MeshcatVisualizer(robot)
while True:
    viz.display(q)
