#!.venv/bin/python3

import time
import numpy as np
import casadi as ca
import pinocchio as pin
from meshcat_viewer_wrapper import MeshcatVisualizer
from pinocchio.robot_wrapper import RobotWrapper

with open('examples/visualization/solution_data/metadata.csv', 'r') as file:
    lines = file.readlines()

# Parse the location
location = lines[0].strip().split(': ')[1].replace("../", "", 2)

# Parse the q0 array
q0_str = lines[1].strip().split(': ')[1]
q0 = np.array([float(x) for x in q0_str.split(', ')])


builder = RobotWrapper.BuildFromURDF
robot = builder(
    location,
    ["."],
    pin.JointModelFreeFlyer(),
)

robot.q0 = q0

# The pinocchio model is what we are really interested by.
model = robot.model
data = model.createData()

ee_names = ["FL_foot", "RL_foot", "FR_foot", "RR_foot"]
start_colors = [(252, 208, 50), (94,75,50), (166,144,50), (177,182,50)]
end_colors = [(252, 208, 255), (94,75,255), (166,144,255), (177,182,255)]


ee_ids = [model.getFrameId(ee_name) for ee_name in ee_names]

new_times = np.genfromtxt('examples/visualization/solution_data/sol_times.csv', delimiter=',')
new_sol = np.genfromtxt('examples/visualization/solution_data/sol_states.csv', delimiter=',')
new_times = np.reshape(np.diff(new_times), (new_times.shape[0] - 1, 1))

# Calculate the number of time steps
num_steps = np.size(new_times, 0)

viz = MeshcatVisualizer(robot)
viz.display(robot.q0)

for k in range(num_steps):
    pin.framesForwardKinematics(model, data, new_sol[:,k])
    ee_positions = {}
    for ee_name in ee_names:
        ee_positions[ee_name] = data.oMf[model.getFrameId(ee_name)].translation

    # Plot the positions of the end effectors
    for i, (ee_name, position) in enumerate(ee_positions.items()):
        # Calculate the current color using linear interpolation
        color = [np.interp(k, [0, num_steps-1], [start_colors[i][j], end_colors[i][j]]) for j in range(3)]
        color = [c/255 for c in color]  # Normalize color

        sphere_name = f"{ee_name}_{k}"  # Append time step to sphere name
        viz.addSphere(sphere_name, radius=0.005, color=color)
        viz.applyConfiguration(sphere_name, list(position) + [0, 0, 0, 1])

def display_scene(q: np.ndarray, dt=1e-1):
    pin.framesForwardKinematics(model, data, q)

    viz.display(q)
    time.sleep(dt)

def display_traj(qs: np.ndarray, dts: np.ndarray):
    for k in range(np.size(dts, 0)):
        display_scene(qs[:, k], dts[k, 0])

while True:
    display_traj(new_sol, new_times)
    time.sleep(1)