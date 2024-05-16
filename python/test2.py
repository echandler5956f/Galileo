import casadi as ca
import numpy as np
import matplotlib.pyplot as plt

# Degree of interpolating polynomial
d = 3

# Get collocation points
tau_root = np.append(0, ca.collocation_points(d, 'radau'))

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
f = ca.Function('f', [x, u], [xdot, L], ['x', 'p'], ['ode', 'quad'])

# Control discretization
N = 1000 # number of control intervals
h = T/N

# Variables for a single collocation interval
Xc = [ca.SX.sym('Xc_'+str(j), 2) for j in range(d)]
X0 = ca.SX.sym('X0', 2)
P = ca.SX.sym('P')
L = ca.SX.sym('L')

# State at the end of the collocation interval
Xf = D[0]*X0
for j in range(d): Xf = Xf + D[j+1]*Xc[j]

# Nonlinear equations and cost contribution for collocation interval
eq = []
Qf = 0
for j in range(d):
    # Expression for the state derivative at the collocation point
    xp = C[0,j+1]*X0
    for r in range(d): xp = xp + C[r+1,j+1]*Xc[r]

    # Evaluate the function at the collocation points
    fj = f(x=Xc[j], p=P)

    # Append collocation equations
    eq.append(h*fj['ode'] - xp)

    # Add cost contribution
    Qf = Qf + B[j+1]*fj['quad']*h

# Implicit discrete-time dynamics
# F = ca.Function('F', [ca.vertcat(*Xc), X0, P], [ca.vertcat(*eq), Xf, Qf],\
#                      ['xc', 'x0', 'p'], ['eq', 'xf', 'qf'])
Feq = ca.Function('feq', [ca.vertcat(*Xc), X0, P], [ca.vertcat(*eq)],\
                     ['xc', 'x0', 'p'], ['eq'])
Fxf = ca.Function('feq', [ca.vertcat(*Xc), X0, P], [Xf],\
                     ['xc', 'x0', 'p'], ['xf'])
Fxq = ca.Function('fxq', [L, ca.vertcat(*Xc), X0, P], [L + Qf])


opti = ca.Opti()
J = 0
var_xs = [opti.variable(2) for k in range(N + 1)]
var_us = [opti.variable(1) for k in range(N)]
var_xcs = [opti.variable(2 * d) for k in range(N)]

Xk = var_xs[0]
opti.subject_to(Xk == [0, 1])
opti.set_initial(Xk, [0, 1])

feqmap = Feq.map(N, "openmp")
fxfmap = Fxf.map(N, "openmp")
fxqfold = Fxq.fold(N)

for k in range(N):
    Uk = var_us[k]
    opti.bounded(-1, Uk, 1)
    Xc = var_xcs[k]
    opti.bounded([-0.25, -np.inf] * d, var_xcs[k], [np.inf, np.inf] * d)
    Xk = var_xs[k + 1]
    opti.bounded([-0.25, -np.inf], Xk, [np.inf, np.inf])

xs = var_xs[0]
us = var_us[0]
xcs = var_xcs[0]
for k in range(1, N):
    xs = ca.horzcat(xs, var_xs[k])
    us = ca.horzcat(us, var_us[k])
    xcs = ca.horzcat(xcs, var_xcs[k])

opti.subject_to(feqmap(xcs, xs, us) == 0)
opti.subject_to(fxfmap(xcs, xs, us) - ca.horzcat(xs[:,1:], var_xs[N]) == 0)
J = fxqfold(J, xcs, xs, us)

opti.minimize(J)
opti.solver("ipopt",  {"expand": True}, {"max_iter": 200, "linear_solver": "ma97"})
# opti.to_function("opti", {var_xs, var_us, var_xcs, opti.lam_g}, {var_xs, var_us, var_xcs})
sol = opti.solve()

# Solve the NLP
x_opt = np.transpose(np.array([opti.value(var_x) for var_x in var_xs]))
u_opt = np.transpose(np.array([opti.value(var_u) for var_u in var_us]))

# Plot the result
tgrid = np.linspace(0, T, N+1)
plt.figure(1)
plt.clf()
plt.plot(tgrid, x_opt[0, :], '--')
plt.plot(tgrid, x_opt[1, :], '-')
plt.step(tgrid, np.append(np.nan, u_opt), '-.')
plt.xlabel('t')
plt.legend(['x1','x2','u'])
plt.grid()
plt.show()

# # # Create a s
# # ll times: [ 0.          1.05662433  3.94337567  5.          6.05662433  8.94337567
# #  10.        ]
# # 7
# all_times = np.zeros((N * (d+1) + 1,1))
# all_times[N * (d+1)] = 10
# for k in range(N):
#     for j in range(d + 1):
#         col_time = tau_root[j] * h + k * h
#         all_times[k * (d + 1) + j] = col_time
# print("All times:", all_times)
# print(all_times.size)