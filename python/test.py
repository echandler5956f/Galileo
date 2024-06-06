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
    

def quat_to_ZYX_euler(quat: np.ndarray) -> np.ndarray:
    qw = quat[3]
    qx = quat[0]
    qy = quat[1]
    qz = quat[2]

    euler = np.zeros((3, 1))
    euler[0] = math.atan2(2 * (qw * qx + qy * qz), 1 - 2 * (qx ** 2 + qy ** 2))
    euler[1] = math.asin(2 * (qw * qy - qz * qx))
    euler[2] = math.atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy ** 2 + qz ** 2))

    return euler


def euler_ZYX_to_quat(eul: np.ndarray) -> np.ndarray:
    cy = math.cos(eul[2] * 0.5)
    sy = math.sin(eul[2] * 0.5)
    cp = math.cos(eul[1] * 0.5)
    sp = math.sin(eul[1] * 0.5)
    cr = math.cos(eul[0] * 0.5)
    sr = math.sin(eul[0] * 0.5)

    qw = cr * cp * cy + sr * sp * sy
    qx = sr * cp * cy - cr * sp * sy
    qy = cr * sp * cy + sr * cp * sy
    qz = cr * cp * sy - sr * sp * cy

    return np.array([[qx], [qy], [qz], [qw]])
    

def eul_int(x: np.ndarray, dx: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    pos = x[:3]
    quat = x[3:]

    eul = quat_to_ZYX_euler(quat)

    vel = dx[:3] * dt
    omega = dx[3:] * dt

    eul_new = eul + angular_vel_to_ZYX_euler_derivative_taylor(omega, eul)
    print(eul_new)

    rodrigues_term = rodrigues_simplified(omega) @ vel
    pos_new = pos + euler_ZYX_to_rot(eul_new) @ rodrigues_term
    quat_new = euler_ZYX_to_quat(eul_new)

    quat_new = quat_new / np.sqrt(
        quat_new[0] ** 2 + quat_new[1] ** 2 + quat_new[2] ** 2 + quat_new[3] ** 2 + 1e-16
    )

    return np.concatenate((pos_new, quat_new), axis=0)
    

def angular_vel_to_ZYX_euler_derivative(omega: np.ndarray, eul: np.ndarray) -> np.ndarray:
    B = np.zeros((3, 3))
    B[0, 1] = np.sin(eul[0])/np.cos(eul[1])
    B[0, 2] = np.cos(eul[0])/np.cos(eul[1])
    B[1, 1] = np.cos(eul[0])
    B[1, 2] = -np.sin(eul[0])
    B[2, 0] = 1
    B[2, 1] = np.sin(eul[0])*np.tan(eul[1])
    B[2, 2] = np.cos(eul[0])*np.tan(eul[1])

    return B @ omega


def angular_vel_to_ZYX_euler_derivative_taylor(omega: np.ndarray, eul: np.ndarray) -> np.ndarray:
    B = np.zeros((3, 3))
    
    # Assuming small angles, we can approximate sin(x) as x and cos(x) as 1
    B[0, 1] = eul[0] / (1 + eul[1]**2 / 2)  # sin(eul[0])/cos(eul[1]) approximated
    B[0, 2] = 1 / (1 + eul[1]**2 / 2)       # cos(eul[0])/cos(eul[1]) approximated
    B[1, 1] = 1                             # cos(eul[0]) approximated
    B[1, 2] = -eul[0]                       # -sin(eul[0]) approximated
    B[2, 0] = 1
    B[2, 1] = eul[0] * eul[1]               # sin(eul[0])*tan(eul[1]) approximated
    B[2, 2] = 1 * eul[1]                    # cos(eul[0])*tan(eul[1]) approximated

    return B @ omega


def euler_ZYX_to_rot(eul: np.ndarray) -> np.ndarray:
    rot = np.zeros((3, 3))
    rot[0, 0] = np.cos(eul[1]) * np.cos(eul[2])
    rot[1, 0] = np.cos(eul[1]) * np.sin(eul[2])
    rot[2, 0] = -np.sin(eul[1])
    rot[0, 1] = np.sin(eul[0]) * np.sin(eul[1]) * np.cos(eul[2]) - np.cos(eul[0]) * np.sin(eul[2])
    rot[1, 1] = np.sin(eul[0]) * np.sin(eul[1]) * np.sin(eul[2]) + np.cos(eul[0]) * np.cos(eul[2])
    rot[2, 1] = np.sin(eul[0]) * np.cos(eul[1])
    rot[0, 2] = np.cos(eul[0]) * np.sin(eul[1]) * np.cos(eul[2]) + np.sin(eul[0]) * np.sin(eul[2])
    rot[1, 2] = np.cos(eul[0]) * np.sin(eul[1]) * np.sin(eul[2]) - np.sin(eul[0]) * np.cos(eul[2])
    rot[2, 2] = np.cos(eul[0]) * np.cos(eul[1])

    return rot
    

def rot_to_quat(R: np.ndarray) -> np.ndarray:
    q = np.zeros((4, 1))
    tr = R[0, 0] + R[1, 1] + R[2, 2]
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2
        q[3] = 0.25 * S
        q[0] = (R[2, 1] - R[1, 2]) / S
        q[1] = (R[0, 2] - R[2, 0]) / S
        q[2] = (R[1, 0] - R[0, 1]) / S
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        S = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2
        q[3] = (R[2, 1] - R[1, 2]) / S
        q[0] = 0.25 * S
        q[1] = (R[0, 1] + R[1, 0]) / S
        q[2] = (R[0, 2] + R[2, 0]) / S
    elif R[1, 1] > R[2, 2]:
        S = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2
        q[3] = (R[0, 2] - R[2, 0]) / S
        q[0] = (R[0, 1] + R[1, 0]) / S
        q[1] = 0.25 * S
        q[2] = (R[1, 2] + R[2, 1]) / S
    else:
        S = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2
        q[3] = (R[1, 0] - R[0, 1]) / S
        q[0] = (R[0, 2] + R[2, 0]) / S
        q[1] = (R[1, 2] + R[2, 1]) / S
        q[2] = 0.25 * S
    return q


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


def quat_to_rot_alt(quat: np.ndarray) -> np.ndarray:
    R_quat = np.zeros((3, 3))

    a = quat[0]
    b = quat[1]
    c = quat[2]
    d = quat[3]

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


def quat_lie_group_int(qb: np.ndarray, vb: np.ndarray, dt: float) -> np.ndarray:
    pos = np.reshape(qb[0:3], (3, 1))
    quat = np.reshape(qb[3:7], (4, 1))

    vel = np.reshape(vb[0:3] * dt, (3, 1))
    omega = np.reshape(vb[3:6] * dt, (3, 1))
    omega_quat = np.concatenate((omega / 2, np.zeros((1, 1))), axis=0)
    exp_omega_quat = quat_exp(omega_quat)
    quat_new = quat_mult(quat, exp_omega_quat)
    quat_new = quat_new / np.sqrt(
        quat_new[0] ** 2 + quat_new[1] ** 2 + quat_new[2] ** 2 + quat_new[3] ** 2 + 1e-16
    )
    pos_new = pos + quat_apply(quat, rodrigues(omega) @ vel)
    return np.concatenate((pos_new, quat_new), axis=0)


def quat_exp(quat: np.ndarray) -> np.ndarray:
    v = quat[0:3]
    v_norm = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + 1e-16)
    a = quat[3]

    exp_a = np.exp(a)
    cos_norm_v = np.cos(v_norm)
    sin_norm_v = np.sin(v_norm)

    quat_new = np.concatenate(
        (exp_a * v * sin_norm_v / v_norm, np.reshape(exp_a * cos_norm_v, (1, 1))),
        axis=0,
    )

    return quat_new / np.sqrt(
        quat_new[0] ** 2 + quat_new[1] ** 2 + quat_new[2] ** 2 + quat_new[3] ** 2 + 1e-16
    )


def quat_apply(quat: np.ndarray, vec3: np.ndarray) -> np.ndarray:
    imag = quat[0:3]
    real = quat[3]
    return vec3 + 2 * np.cross(imag, np.cross(imag, vec3, axis=0) + real * vec3, axis=0)


def quat_mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    result = np.zeros((4, 1))
    real_1 = q1[3]
    imag_1 = q1[0:3].reshape(-1)

    real_2 = q2[3]
    imag_2 = q2[0:3].reshape(-1)

    real = real_1 * real_2 - np.dot(imag_1, imag_2)
    imag = real_1 * imag_2 + real_2 * imag_1 + np.cross(imag_1, imag_2, axis=0)

    result[0:3] = imag.reshape((3, 1))
    result[3] = real

    return result


def quat_inv(q: np.ndarray) -> np.ndarray:
    q_norm = np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2 + 1e-16)
    return np.concatenate((-q[0:3] / q_norm, q[3] / q_norm), axis=0)


def quat_distance(quat1: np.ndarray, quat2: np.ndarray) -> float:
    real_1 = quat1[3]
    imag_1 = quat1[0:3]

    real_2 = quat2[3]
    imag_2 = quat2[0:3]

    return real_1 * imag_2 - real_2 * imag_1 + np.cross(imag_1, imag_2, axis=0)


def quat_distance_alt(quat1: np.ndarray, quat2: np.ndarray) -> float:
    real_1 = quat1[0]
    imag_1 = quat1[1:]

    real_2 = quat2[0]
    imag_2 = quat2[1:]

    return real_1 * imag_2 - real_2 * imag_1 + np.cross(imag_1, imag_2, axis=0)


def rodrigues(omega: np.ndarray) -> np.ndarray:
    theta = np.sqrt(omega[0] ** 2 + omega[1] ** 2 + omega[2] ** 2 + 1e-16)
    return (
        np.eye(3)
        + ((1 - np.cos(theta)) / (theta * theta)) * skew(omega)
        + ((theta - np.sin(theta)) / (theta * theta * theta))
        * (skew(omega) @ skew(omega))
    )

# def taylor_approximation(omega):
#     theta = np.sqrt(omega[0] ** 2 + omega[1] ** 2 + omega[2] ** 2 + 1e-16)
    
#     # Taylor series expansion for small theta
#     cos_theta_approx = 1 - theta**2 / 2
#     sin_theta_approx = theta - theta**3 / 6
    
#     return (
#         np.eye(3)
#         + ((1 - cos_theta_approx) / (theta * theta)) * skew(omega)
#         + ((theta - sin_theta_approx) / (theta * theta * theta))
#         * (skew(omega) @ skew(omega))
#     )

def rodrigues_simplified(omega: np.ndarray) -> np.ndarray:
    theta = np.sqrt(omega[0] ** 2 + omega[1] ** 2 + omega[2] ** 2 + 1e-16)

    # Taylor series expansion for small theta
    cos_theta_approx = 1 - theta**2 / 2
    sin_theta_approx = theta - theta**3 / 6

    return (
        np.eye(3)
        + ((1 - cos_theta_approx) / (theta * theta)) * skew(omega)
        + ((theta - sin_theta_approx) / (theta * theta * theta))
        * (skew(omega) @ skew(omega))
    )


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


def normalize_rotation_matrix(R: np.ndarray) -> np.ndarray:
    # Gram-Schmidt process
    U, _, Vt = np.linalg.svd(R)
    return U @ Vt


def alt_int(x: np.ndarray, dx: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    pos = x[0:3]
    quat = x[3:7]

    vel = dx[0:3] * dt
    omega = dx[3:6] * dt

    rot_quat = quat_to_rot(quat)
    exp_omega_rot = rot_exp_taylor_series(skew(omega), 8)
    rot_new = rot_quat @ exp_omega_rot
    rot_new = normalize_rotation_matrix(rot_new)  # normalize the rotation matrix
    t = pos + rot_quat @ rodrigues(omega) @ vel
    quat_new = rot_to_quat(rot_new)
    return np.concatenate((t, quat_new), axis=0)


# def main():
#     # Define the initial state
#     x = np.reshape(np.array([0, 0, 0, 0, 0, 0, 1]), (7, 1))
#     dx = np.reshape(np.array([0, 0, 0, 0, -2, 0.1]), (6, 1))

#     # Define the time step
#     dt = 3

#     # Integrate the state using the Lie group integration
#     x_new = quat_lie_group_int(x, dx, dt)
#     # Convert the quaternion to a rotation matrix
#     R = quat_to_rot(x_new[3:7])
#     print(x_new[0:3])
#     print(R)

#     # Integrate the state using the alternative integration
#     t, rot_new = alt_int(x, dx, dt)
#     print(t)
#     print(rot_new)


# if __name__ == "__main__":
#     main()


import numpy as np
import matplotlib.pyplot as plt
from pytransform3d.rotations import (
    quaternion_integrate,
    matrix_from_quaternion,
    plot_basis,
)
from pytransform3d.rotations import q_id
from pytransform3d.batch_rotations import quaternion_slerp_batch
from pytransform3d.trajectories import plot_trajectory

N = 50
dt = 0.025

x0 = np.reshape(np.array([0, 0, 0, 0, 0, 0, 1]), (7, 1))
# velocities = np.empty((N, 6))
# velocities[:, :] = np.array([0, 0, 0, 0, 0, 0])

np.random.seed(35)
velocities = 5 * np.random.random((N, 6)) - 2.5
# velocities[:, :] = np.array([0.1, 0.1, 0.1, 0, 0, (np.pi / 2) / (N * dt)])
# randvelocity = 0.2 * np.random.random((1, 6)) - 0.1
# velocities[:,:] = randvelocity
# velocities[N-1, :] = np.array([0, 0, 0, 0, 0, 0])

# for t in range(N-1):
#     velocities[t+1, :] = velocities[t, :] + np.array([0, 0, 0, 0, 0, 0.1])

# print(velocities)

Q1s = quaternion_integrate(velocities[:, 3:], dt=dt)
P1s = np.empty((N, 3))

x = x0.copy()
P1s[0, :] = x[:3].flatten()
for t in range(N - 1):
    vel = velocities[t, :3] * dt
    omega = velocities[t, 3:] * dt
    exp_rot = quat_to_rot_alt(Q1s[t])
    Vt = rodrigues(omega) @ vel
    M_exp = np.eye(4)
    M_exp[:3, :3] = exp_rot
    M_exp[:3, 3] = Vt.flatten()

    M0 = np.eye(4)
    M0[:3, :3] = exp_rot
    M0[:3, 3] = P1s[t, :]

    M1 = M0 @ M_exp
    P1s[t + 1] = M1[:3, 3]


x = x0.copy()
Q2s = np.empty((N, 4))
P2s = np.empty((N, 3))

tmp_x = x0.copy()
tmp_x[3] = x[6]
tmp_x[4] = x[3]
tmp_x[5] = x[4]
tmp_x[6] = x[5]
Q2s[0] = tmp_x[3:].flatten()
P2s[0] = tmp_x[:3].flatten()
for t in range(1, N):
    dx = velocities[t-1].reshape((6, 1))
    x = quat_lie_group_int(x, dx, dt)
    tmp_x = x.copy()
    tmp_x[3] = x[6]
    tmp_x[4] = x[3]
    tmp_x[5] = x[4]
    tmp_x[6] = x[5]
    Q2s[t] = tmp_x[3:].flatten()
    P2s[t] = tmp_x[:3].flatten()


x = x0.copy()
Q3s = np.empty((N, 4))
P3s = np.empty((N, 3))

tmp_x = x0.copy()
tmp_x[3] = x[6]
tmp_x[4] = x[3]
tmp_x[5] = x[4]
tmp_x[6] = x[5]
Q3s[0] = tmp_x[3:].flatten()
P3s[0] = tmp_x[:3].flatten()
for t in range(1, N):
    dx = velocities[t-1].reshape((6, 1))
    x = alt_int(x, dx, dt)
    tmp_x = x.copy()
    tmp_x[3] = x[6]
    tmp_x[4] = x[3]
    tmp_x[5] = x[4]
    tmp_x[6] = x[5]
    Q3s[t] = tmp_x[3:].flatten()
    P3s[t] = tmp_x[:3].flatten()


x = x0.copy()
Q4s = np.empty((N, 4))
P4s = np.empty((N, 3))

tmp_x = x0.copy()
tmp_x[3] = x[6]
tmp_x[4] = x[3]
tmp_x[5] = x[4]
tmp_x[6] = x[5]
Q4s[0] = tmp_x[3:].flatten()
P4s[0] = tmp_x[:3].flatten()
for t in range(1, N):
    dx = velocities[t-1].reshape((6, 1))
    x = eul_int(x, dx, dt)
    tmp_x = x.copy()
    tmp_x[3] = x[6]
    tmp_x[4] = x[3]
    tmp_x[5] = x[4]
    tmp_x[6] = x[5]
    Q4s[t] = tmp_x[3:].flatten()
    P4s[t] = tmp_x[:3].flatten()

print("--------------------------------------------------\n")
# print(Q1s)
# print(Q2s)
# print(Q3s)

print("--------------------------------------------------\n")

# print(P1s)
# print(P2s)
# print(P3s)
# print("--------------------------------------------------\n")

# ax = None
# for t in range(N):
#     R = matrix_from_quaternion(Q1s[t])
#     P = P1s[t]
#     ax = plot_basis(ax=ax, s=0.15, R=R, p=P)
# plt.show()

# ax = None
# for t in range(N):
#     R = matrix_from_quaternion(Q2s[t])
#     P = P2s[t]
#     ax = plot_basis(ax=ax, s=0.15, R=R, p=P)
# plt.show()

error_p = 0
error_q = 0
for t in range(N):
    error_p = (
        error_p
        + (P1s[t, 0] - P4s[t, 0]) ** 2 #+ (P2s[t, 0] - P3s[t, 0]) ** 2
        + (P1s[t, 1] - P4s[t, 1]) ** 2 #+ (P2s[t, 1] - P3s[t, 1]) ** 2
        + (P1s[t, 2] - P4s[t, 2]) ** 2 #+ (P2s[t, 2] - P3s[t, 2]) ** 2
    )
    quat_dist1 = quat_distance_alt(Q1s[t], Q4s[t])
    quat_dist2 = quat_distance_alt(Q1s[t], Q4s[t])
    error_q = error_q + quat_dist1[0] ** 2 + quat_dist1[1] ** 2 + quat_dist1[2] ** 2 
    # error_q = error_q + quat_dist2[0] ** 2 + quat_dist2[1] ** 2 + quat_dist2[2] ** 2

error_p = np.sqrt(error_p)
error_q = np.sqrt(error_q)

print(error_p)
print(error_q)
