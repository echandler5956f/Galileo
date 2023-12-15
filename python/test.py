import casadi as ca
import numpy as np
import matplotlib.pyplot as plt

# Degree of interpolating polynomial
d = 3

# Get collocation points
tau_root = np.append(0, ca.collocation_points(d, 'legendre'))

# Coefficients of the collocation equation
C = np.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = np.zeros(d+1)

# Coefficients of the quadrature function
B = np.zeros(d+1)

# Construct polynomial basis
for j in range(d+1):
    # Construct Lagrange polynomials to get the polynomial basis at the collocation point
    p = np.poly1d([1])
    for r in range(d+1):
        if r != j:
            p *= np.poly1d([1, -tau_root[r]]) / (tau_root[j]-tau_root[r])

    # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    D[j] = p(1.0)

    # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    pder = np.polyder(p)
    for r in range(d+1):
        C[j,r] = pder(tau_root[r])

    # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
    pint = np.polyint(p)
    B[j] = pint(1.0)

# Time horizon
T = 10.

# Declare model variables
x1 = ca.SX.sym('x1')
x2 = ca.SX.sym('x2')
x = ca.vertcat(x1, x2)
u = ca.SX.sym('u')

# Model equations
xdot = ca.vertcat((1-x2**2)*x1 - x2 + u, x1)

# Objective term
L = x1**2 + x2**2 + u**2

# Continuous time dynamics
f = ca.Function('f', [x, u], [xdot, L], ['x', 'u'], ['xdot', 'L'])

# Control discretization
N = 20  # Number of control intervals
h = T/N

# Start with an empty NLP
w = []
w0 = []
lbw = []
ubw = []
J = 0
g = []
lbg = []
ubg = []

# For plotting x and u given w
x_plot = []
x_all = []
u_plot = []

# "Lift" initial conditions
Xk = ca.MX.sym('X0', 2)
w.append(Xk)
lbw.append([0, 1])
ubw.append([0, 1])
w0.append([0, 1])
x_plot.append(Xk)
x_all.append(Xk)

# Formulate the NLP
for k in range(N):
    # New NLP variable for the control
    Uk = ca.MX.sym('U_' + str(k))
    w.append(Uk)
    lbw.append([-1])
    ubw.append([1])
    w0.append([0])
    u_plot.append(Uk)

    # State at Chebyshev collocation points
    Xc = []
    for j in range(d):
        Xkj = ca.MX.sym('X_' + str(k) + '_' + str(j), 2)
        x_all.append(Xkj)
        Xc.append(Xkj)
        w.append(Xkj)
        lbw.append([-0.25, -np.inf])
        ubw.append([np.inf,  np.inf])
        w0.append([0, 0])

    # Loop over Chebyshev collocation points
    Xk_end = D[0]*Xk
    for j in range(1, d+1):
        # Expression for the state derivative at the Chebyshev collocation point
        xp = C[0, j]*Xk
        for r in range(d):
            xp = xp + C[r+1, j]*Xc[r]

        # Append collocation equations
        fj, qj = f(Xc[j-1], Uk)
        g.append(h*fj - xp)
        lbg.append([0, 0])
        ubg.append([0, 0])

        # Add contribution to the end state
        Xk_end = Xk_end + D[j]*Xc[j-1]

        # Add contribution to the quadrature function
        J = J + B[j]*qj*h

    # New NLP variable for state at the end of the interval
    Xk = ca.MX.sym('X_' + str(k+1), 2)
    w.append(Xk)
    lbw.append([-0.25, -np.inf])
    ubw.append([np.inf,  np.inf])
    w0.append([0, 0])
    x_plot.append(Xk)
    if k < N-1:
        x_all.append(Xk)

    # Add equality constraint
    g.append(Xk_end - Xk)
    lbg.append([0, 0])
    ubg.append([0, 0])

# Concatenate vectors
w = ca.vertcat(*w)
g = ca.vertcat(*g)
x_plot = ca.horzcat(*x_plot)
x_all = ca.horzcat(*x_all)
u_plot = ca.horzcat(*u_plot)
w0 = np.concatenate(w0)
lbw = np.concatenate(lbw)
ubw = np.concatenate(ubw)
lbg = np.concatenate(lbg)
ubg = np.concatenate(ubg)

# Create an NLP solver
prob = {'f': J, 'x': w, 'g': g}
solver = ca.nlpsol('solver', 'ipopt', prob)

# Function to get x and u trajectories from w
trajectories = ca.Function('trajectories', [w], [x_plot, u_plot], ['w'], ['x', 'u'])
trajectoriesAll = ca.Function('trajectoriesAll', [w], [x_all], ['w'], ['x_all'])

# Solve the NLP
sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
x_opt, u_opt = trajectories(sol['x'])
xall_opt = trajectoriesAll(sol['x'])
x_opt = x_opt.full()  # to numpy array
xall_opt = xall_opt.full()
u_opt = u_opt.full()  # to numpy array

xall_opt1 = np.reshape(xall_opt[0], (N, d + 1))
xall_opt2 = np.reshape(xall_opt[1], (N, d + 1))

from scipy.interpolate import BarycentricInterpolator

# Plot the result
tgrid = np.linspace(0, T, N + 1)

# Create barycentric interpolators for x1 and x2 for each segment
interpolators_x1 = []
interpolators_x2 = []

for k in range(N):
    t_segment = np.linspace(k * h, (k + 1) * h, d + 1)  # Duration of each segment
    x1_segment = xall_opt1[k]
    x2_segment = xall_opt2[k]

    x1_interp = BarycentricInterpolator(t_segment, x1_segment)
    x2_interp = BarycentricInterpolator(t_segment, x2_segment)

    interpolators_x1.append(x1_interp)
    interpolators_x2.append(x2_interp)

# Create a dense grid for interpolation over the entire trajectory
t_dense = np.linspace(0, T, 2 * N * (d + 1))

# Interpolate the values for the entire trajectory
x_interp_values_x1 = np.concatenate([np.hstack((x1_optk, interp(t_segment))) for x1_optk, interp, t_segment in zip(x_opt[0], interpolators_x1, np.array_split(t_dense, N))])
x_interp_values_x2 = np.concatenate([np.hstack((x2_optk, interp(t_segment))) for x2_optk, interp, t_segment in zip(x_opt[1], interpolators_x2, np.array_split(t_dense, N))])
t_dense = np.linspace(0, T, 2 * N * (d + 1) + N)

# Plot the interpolated values
plt.figure(1)
plt.clf()
plt.step(tgrid, np.append(np.nan, u_opt[0]), '-.')
plt.plot(tgrid, x_opt[0], '.')
plt.plot(t_dense, x_interp_values_x1, '-', label='x1')
# plt.plot(t_dense, x_interp_values_x2, '--', label='x2')
plt.xlabel('t')
plt.legend()
plt.grid()
plt.show()

# # Create a dense grid for interpolation over the entire trajectory
# t_dense = np.linspace(0, T, 2 * N * (d + 1))

# # Interpolate the values for the entire trajectory
# x_interp_values_x1 = np.concatenate([interp(t_segment) for interp, t_segment in zip(interpolators_x1, np.array_split(t_dense, N))])
# x_interp_values_x2 = np.concatenate([interp(t_segment) for interp, t_segment in zip(interpolators_x2, np.array_split(t_dense, N))])

# # Plot the interpolated values
# plt.figure(1)
# plt.clf()
# plt.step(tgrid, np.append(np.nan, u_opt[0]), '-.')
# plt.plot(tgrid, x_opt[0], '.')
# plt.plot(tgrid, x_opt[1], '.')
# plt.plot(t_dense, x_interp_values_x1, 'o-', label='x1 (interpolated)')
# plt.plot(t_dense, x_interp_values_x2, 'o-', label='x2 (interpolated)')


# # Plot the knot points (boundaries of each segment)
# for k in range(N):
#     t_knot_start = k * T / N
#     t_knot_end = (k + 1) * T / N
#     x1_knot_start = x_opt[0][k]
#     x1_knot_end = x_opt[0][k + 1]
    
#     # plt.plot([t_knot_start, t_knot_start], [x1_knot_start, x_interp_values_x1[(k+1)*(d+1)-1]], 'ro', label=f'x1 (knot start) segment {k + 1}')
#     # plt.plot([t_knot_end, t_knot_end], [x1_knot_end, x_interp_values_x1[(k+1)*(d+1)]], 'ro', label=f'x1 (knot end) segment {k + 1}')

# plt.xlabel('t')
# plt.legend()
# plt.grid()
# plt.show()