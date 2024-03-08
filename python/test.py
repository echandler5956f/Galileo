# import casadi as ca
# import numpy as np
# import matplotlib.pyplot as plt


# # Rosenbrock
# # nx = 2
# # nu = 1

# # dx = 9
# # du = 1
# # N = 2
# # phases = 2

# # num_variables = (nx*(dx*N + N + 1) + nu*(du*N + N + 1)) * phases
# # print(num_variables)

# # phase_continuity_constraints = nx * phases
# # collocation_constraints = nx * dx * N * phases
# # final_state_constraint = nx * N * phases
# # final_input_constraint = nu * N * phases

# # num_equality_constraints = phase_continuity_constraints + collocation_constraints + final_state_constraint + final_input_constraint
# # print(num_equality_constraints)

# # Go1 Quadruped
# nx = 49
# ndx = 48
# nu = 24

# num_feet = 4

# dx = 3
# du = 2
# N = 2
# feet_con = [num_feet, 2, num_feet]
# phases = len(feet_con)

# num_variables = (ndx*(dx*N + N + 1) + nu*(du*N + N + 1)) * phases
# print(num_variables)

# phase_continuity_constraints = ndx * (phases - 1) + nx
# collocation_constraints = ndx * dx * N * phases
# final_state_constraint = ndx * N * phases
# final_input_constraint = nu * N * phases
# velocity_constraint = 0
# force_constraint = 0
# for i in range(phases):
#     velocity_constraint += 3 * feet_con[i] * dx * N
#     force_constraint += 3 * (num_feet - feet_con[i]) * dx * N

# num_equality_constraints = phase_continuity_constraints + collocation_constraints + final_state_constraint + final_input_constraint + velocity_constraint + force_constraint
# print(num_equality_constraints)

# # # Degree of interpolating polynomial
# # d = 3

# # # Get collocation points
# # tau_root = np.append(0, ca.collocation_points(d, 'legendre'))

# # # Coefficients of the collocation equation
# # C = np.zeros((d+1,d+1))

# # # Coefficients of the continuity equation
# # D = np.zeros(d+1)

# # # Coefficients of the quadrature function
# # B = np.zeros(d+1)

# # # Construct polynomial basis
# # for j in range(d+1):
# #     # Construct Lagrange polynomials to get the polynomial basis at the collocation point
# #     p = np.poly1d([1])
# #     for r in range(d+1):
# #         if r != j:
# #             p *= np.poly1d([1, -tau_root[r]]) / (tau_root[j]-tau_root[r])

# #     # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
# #     D[j] = p(1.0)

# #     # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
# #     pder = np.polyder(p)
# #     for r in range(d+1):
# #         C[j,r] = pder(tau_root[r])

# #     # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
# #     pint = np.polyint(p)
# #     B[j] = pint(1.0)

# # # Time horizon
# # T = 10.

# # # Declare model variables
# # x1 = ca.SX.sym('x1')
# # x2 = ca.SX.sym('x2')
# # x = ca.vertcat(x1, x2)
# # u = ca.SX.sym('u')

# # # Model equations
# # xdot = ca.vertcat((1-x2**2)*x1 - x2 + u, x1)

# # # Objective term
# # L = x1**2 + x2**2 + u**2

# # # Continuous time dynamics
# # f = ca.Function('f', [x, u], [xdot, L], ['x', 'u'], ['xdot', 'L'])

# # # Control discretization
# # N = 5  # Number of control intervals
# # h = T/N

# # # Start with an empty NLP
# # w = []
# # w0 = []
# # lbw = []
# # ubw = []
# # J = 0
# # g = []
# # lbg = []
# # ubg = []

# # # For plotting x and u given w
# # x_plot = []
# # x_all = []
# # u_plot = []

# # # "Lift" initial conditions
# # Xk = ca.MX.sym('X0', 2)
# # w.append(Xk)
# # lbw.append([0, 1])
# # ubw.append([0, 1])
# # w0.append([0, 1])
# # x_plot.append(Xk)
# # x_all.append(Xk)

# # # Formulate the NLP
# # for k in range(N):
# #     # New NLP variable for the control
# #     Uk = ca.MX.sym('U_' + str(k))
# #     w.append(Uk)
# #     lbw.append([-1])
# #     ubw.append([1])
# #     w0.append([0])
# #     u_plot.append(Uk)

# #     # State at Chebyshev collocation points
# #     Xc = []
# #     for j in range(d):
# #         Xkj = ca.MX.sym('X_' + str(k) + '_' + str(j), 2)
# #         x_all.append(Xkj)
# #         Xc.append(Xkj)
# #         w.append(Xkj)
# #         lbw.append([-0.25, -np.inf])
# #         ubw.append([np.inf,  np.inf])
# #         w0.append([0, 0])

# #     # Loop over Chebyshev collocation points
# #     Xk_end = D[0]*Xk
# #     for j in range(1, d+1):
# #         # Expression for the state derivative at the Chebyshev collocation point
# #         xp = C[0, j]*Xk
# #         for r in range(d):
# #             xp = xp + C[r+1, j]*Xc[r]

# #         # Append collocation equations
# #         fj, qj = f(Xc[j-1], Uk)
# #         g.append(h*fj - xp)
# #         lbg.append([0, 0])
# #         ubg.append([0, 0])

# #         # Add contribution to the end state
# #         Xk_end = Xk_end + D[j]*Xc[j-1]

# #         # Add contribution to the quadrature function
# #         J = J + B[j]*qj*h

# #     # New NLP variable for state at the end of the interval
# #     Xk = ca.MX.sym('X_' + str(k+1), 2)
# #     w.append(Xk)
# #     lbw.append([-0.25, -np.inf])
# #     ubw.append([np.inf,  np.inf])
# #     w0.append([0, 0])
# #     x_plot.append(Xk)
# #     if k < N-1:
# #         x_all.append(Xk)

# #     # Add equality constraint
# #     g.append(Xk_end - Xk)
# #     lbg.append([0, 0])
# #     ubg.append([0, 0])

# # # Concatenate vectors
# # w = ca.vertcat(*w)
# # g = ca.vertcat(*g)
# # x_plot = ca.horzcat(*x_plot)
# # x_all = ca.horzcat(*x_all)
# # u_plot = ca.horzcat(*u_plot)
# # w0 = np.concatenate(w0)
# # lbw = np.concatenate(lbw)
# # ubw = np.concatenate(ubw)
# # lbg = np.concatenate(lbg)
# # ubg = np.concatenate(ubg)

# # # Create an NLP solver
# # prob = {'f': J, 'x': w, 'g': g}
# # solver = ca.nlpsol('solver', 'ipopt', prob)

# # # Function to get x and u trajectories from w
# # trajectories = ca.Function('trajectories', [w], [x_plot, u_plot], ['w'], ['x', 'u'])
# # trajectoriesAll = ca.Function('trajectoriesAll', [w], [x_all], ['w'], ['x_all'])

# # # Solve the NLP
# # sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
# # x_opt, u_opt = trajectories(sol['x'])
# # xall_opt = trajectoriesAll(sol['x'])
# # x_opt = x_opt.full()  # to numpy array
# # xall_opt = xall_opt.full()
# # u_opt = u_opt.full()  # to numpy array

# # xall_opt1 = np.reshape(xall_opt[0], (N, d + 1))
# # xall_opt2 = np.reshape(xall_opt[1], (N, d + 1))

# # from scipy.interpolate import BarycentricInterpolator

# # # Plot the result
# # tgrid = np.linspace(0, T, N + 1)

# # # Create barycentric interpolators for x1 and x2 for each segment
# # interpolators_x1 = []
# # interpolators_x2 = []

# # for k in range(N):
# #     t_segment = np.linspace(k * h, (k + 1) * h, d + 1)  # Duration of each segment
# #     x1_segment = xall_opt1[k]
# #     x2_segment = xall_opt2[k]

# #     x1_interp = BarycentricInterpolator(t_segment, x1_segment)
# #     x2_interp = BarycentricInterpolator(t_segment, x2_segment)

# #     interpolators_x1.append(x1_interp)
# #     interpolators_x2.append(x2_interp)

# # # Create a dense grid for interpolation over the entire trajectory
# # t_dense = np.linspace(0, T, 2 * N * (d + 1))

# # # Interpolate the values for the entire trajectory
# # x_interp_values_x1 = np.concatenate([np.hstack((x1_optk, interp(t_segment))) for x1_optk, interp, t_segment in zip(x_opt[0], interpolators_x1, np.array_split(t_dense, N))])
# # x_interp_values_x2 = np.concatenate([np.hstack((x2_optk, interp(t_segment))) for x2_optk, interp, t_segment in zip(x_opt[1], interpolators_x2, np.array_split(t_dense, N))])
# # t_dense = np.linspace(0, T, 2 * N * (d + 1) + N)

# # # Plot the interpolated values
# # plt.figure(1)
# # plt.clf()
# # plt.step(tgrid, np.append(np.nan, u_opt[0]), '-.', label='u')
# # plt.plot(tgrid, x_opt[0], '.')
# # plt.plot(tgrid, x_opt[1], '.')
# # plt.plot(t_dense, x_interp_values_x1, '-', label='x1')
# # plt.plot(t_dense, x_interp_values_x2, '--', label='x2')
# # plt.xlabel('t')
# # plt.legend()
# # plt.grid()
# # plt.show()

# # # # Create a dense grid for interpolation over the entire trajectory
# # # t_dense = np.linspace(0, T, 2 * N * (d + 1))

# # # # Interpolate the values for the entire trajectory
# # # x_interp_values_x1 = np.concatenate([interp(t_segment) for interp, t_segment in zip(interpolators_x1, np.array_split(t_dense, N))])
# # # x_interp_values_x2 = np.concatenate([interp(t_segment) for interp, t_segment in zip(interpolators_x2, np.array_split(t_dense, N))])

# # # # Plot the interpolated values
# # # plt.figure(1)
# # # plt.clf()
# # # plt.step(tgrid, np.append(np.nan, u_opt[0]), '-.')
# # # plt.plot(tgrid, x_opt[0], '.')
# # # plt.plot(tgrid, x_opt[1], '.')
# # # plt.plot(t_dense, x_interp_values_x1, 'o-', label='x1 (interpolated)')
# # # plt.plot(t_dense, x_interp_values_x2, 'o-', label='x2 (interpolated)')


# # # # Plot the knot points (boundaries of each segment)
# # # for k in range(N):
# # #     t_knot_start = k * T / N
# # #     t_knot_end = (k + 1) * T / N
# # #     x1_knot_start = x_opt[0][k]
# # #     x1_knot_end = x_opt[0][k + 1]

# # #     # plt.plot([t_knot_start, t_knot_start], [x1_knot_start, x_interp_values_x1[(k+1)*(d+1)-1]], 'ro', label=f'x1 (knot start) segment {k + 1}')
# # #     # plt.plot([t_knot_end, t_knot_end], [x1_knot_end, x_interp_values_x1[(k+1)*(d+1)]], 'ro', label=f'x1 (knot end) segment {k + 1}')

# # # plt.xlabel('t')
# # # plt.legend()
# # # plt.grid()
# # # plt.show()


# Import numpy for array operations and math for trigonometric functions
import numpy as np
import math
from typing import Tuple


def skew(x):
    if isinstance(x, np.ndarray) and len(x.shape) >= 2:
        return np.array(
            [[0, -x[2][0], x[1][0]], [x[2][0], 0, -x[0][0]], [-x[1][0], x[0][0], 0]]
        )
    else:
        return np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])


def quat_to_rot(quat: np.ndarray) -> np.ndarray:
    R_quat = np.zeros((3, 3))

    a = quat[3]
    b = quat[0]
    c = quat[1]
    d = quat[2]

    s = 2 / (a * a + b * b + c * c + d * d)
    bs = b * s
    cs = c * s
    ds = d * s
    ab = a * bs
    ac = a * cs
    ad = a * ds
    bb = b * bs
    bc = b * cs
    bd = b * ds
    cc = c * cs
    cd = c * ds
    dd = d * ds

    R_quat[0, 0] = 1 - cc - dd
    R_quat[0, 1] = bc - ad
    R_quat[0, 2] = bd + ac
    R_quat[1, 0] = bc + ad
    R_quat[1, 1] = 1 - bb - dd
    R_quat[1, 2] = cd - ab
    R_quat[2, 0] = bd - ac
    R_quat[2, 1] = cd + ab
    R_quat[2, 2] = 1 - bb - cc
    return R_quat


def SE3_from_Rt(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    T = np.zeros((4, 4))
    T[0:3, 0:3] = R
    T[0:3, 3] = t
    T[3, 3] = 1
    return T


def apply_quat(quat: np.ndarray, vec3: np.ndarray) -> np.ndarray:
    imag = quat[0:3]
    real = quat[3]
    return vec3 + 2 * np.cross(imag, np.cross(imag, vec3) + real * vec3)


def omega_to_unit_quat(omega: np.ndarray) -> np.ndarray:
    quat = np.zeros((4, 1))
    t = np.sqrt(omega[0] ** 2 + omega[1] ** 2 + omega[2] ** 2 + 1e-8)
    s = math.sin(t / 2) / t

    quat[0:3] = omega * s
    quat[3] = math.cos(t / 2)
    return quat


def quat_mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    result = np.zeros((4, 1))
    real_1 = q1[3]
    imag_1 = q1[0:3].reshape(-1)

    real_2 = q2[3]
    imag_2 = q2[0:3].reshape(-1)

    real = real_1 * real_2 - np.dot(imag_1, imag_2)
    imag = real_1 * imag_2 + real_2 * imag_1 + np.cross(imag_1, imag_2)

    result[0:3] = imag.reshape((3, 1))
    result[3] = real

    return result


# Define a function that takes a quaternion as an input and returns its exponential as an output
def quat_exp(q: np.ndarray) -> np.ndarray:
    # Extract the scalar and vector parts of the quaternion
    a = q[3]  # scalar part
    v = q[:3]  # vector part

    # Compute the norm of the vector part
    norm_v = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + 1e-8)

    # Compute the exponential of the scalar part
    exp_a = math.exp(a)

    # Compute the cosine and sine of the norm of the vector part
    cos_norm_v = math.cos(norm_v)
    sin_norm_v = math.sin(norm_v)

    # Compute the exponential of the quaternion using the formula
    exp_q = np.zeros((4, 1))  # initialize an array of zeros
    exp_q[3] = exp_a * cos_norm_v  # scalar part
    exp_q[:3] = exp_a * (v / norm_v) * sin_norm_v  # vector part

    # Normalize the exponential of the quaternion
    exp_q = exp_q / np.sqrt(
        exp_q[0] ** 2 + exp_q[1] ** 2 + exp_q[2] ** 2 + exp_q[3] ** 2 + 1e-8
    )

    # Return the exponential of the quaternion
    return exp_q


def rot_exp_taylor_series(rot: np.ndarray, degree: int) -> np.ndarray:
    # initialize the result as the identity matrix
    result = np.eye(3)
    # initialize the power of the matrix as the matrix itself
    power = rot
    # initialize the factorial as 1
    factorial = 1
    # loop from 1 to d
    for i in range(1, degree + 1):
        # add the power of the matrix divided by the factorial to the result
        result = result + power / factorial
        # multiply the power of the matrix by the matrix
        power = power @ rot
        # increment the factorial
        factorial = factorial * (i + 1)
    # return the result
    return result


def lie_group_int(x: np.ndarray, dx: np.ndarray, dt: float) -> np.ndarray:
    pos = x[0:3]
    quat = x[3:7]

    vel = dx[0:3]
    omega = dx[3:6]

    q_omega = omega_to_unit_quat(omega)
    exp_q = quat_exp(q_omega * dt)
    q_new = quat_mult(exp_q, quat)

    R = quat_to_rot(q_new)
    t = pos + vel * dt
    return np.concatenate((t, q_new), axis=0)


def alt_int(x: np.ndarray, dx: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    pos = x[0:3]
    quat = x[3:7]

    vel = dx[0:3]
    omega = dx[3:6]

    rot_quat = quat_to_rot(quat)
    exp_omega_rot = rot_exp_taylor_series(skew(omega), 8)
    rot_new = rot_quat @ exp_omega_rot
    t = pos + vel * dt
    return (t, rot_new)


def main():
    # Define the initial state
    x = np.reshape(np.array([0, 0, 0, 0, 0, 0, 1]), (7, 1))
    dx = np.reshape(np.array([1, 0, 0, 0, 0, 0]), (6, 1))

    # Define the time step
    dt = 1.0

    # Integrate the state using the Lie group integration
    x_new = lie_group_int(x, dx, dt)
    # Convert the quaternion to a rotation matrix
    R = quat_to_rot(x_new[3:7])
    print(x_new[0:3])
    print(R)

    # Integrate the state using the alternative integration
    t, rot_new = alt_int(x, dx, dt)
    print(t)
    print(rot_new)


if __name__ == "__main__":
    main()
