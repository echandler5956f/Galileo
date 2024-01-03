#!.venv/bin/python3

from frozendict import frozendict
from collections import defaultdict
from typing import List
import time
import numpy as np
import casadi as ca
import pinocchio as pin
from dataclasses import dataclass
from pinocchio import casadi as cpin
from meshcat_viewer_wrapper import MeshcatVisualizer
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation as R
from scipy.spatial.transform import Slerp


# Class for storing end efector information
@dataclass(frozen=True)
class EndEffector:
    # Name of the frame in the URDF
    frame_name: str
    # ID of the frame in pinocchio
    frame_id: int
    # True if the end effector is 6-DOF, False if 3-DOF
    type_6D: bool


# Class for storing phase
@dataclass(frozen=True)
class Phase:
    # Dictionary which maps end effectors to contact flags (True is in contact)
    contacts: frozendict({EndEffector: bool})
    # Name of the phase (used for hashing)
    phase_name: str = None
    # If the phase has a fixed period, this value is set
    fixed_timing: float = None
    # If the phase has a variable period (determined using Casadi), this flag is True
    timing_var: bool = False
    # Number of knot points in this phase
    knot_points: int = 1


# Class for operating on contacts
class ContactSequence:
    def __init__(self):
        self.sequence: List[Phase] = []
        self.cumulative_knots: List[int] = []

    # Add a new phase to the sequence
    # Cumulative_knots[-1] returns the total number of knot points in the sequence
    def add_phase(self, phase: Phase):
        self.sequence.append(phase)

        previous_knots = self.cumulative_knots[-1] if self.cumulative_knots else 0
        self.cumulative_knots.append(previous_knots + phase.knot_points)

    # Get a list of the end effectors (useful for iterating over)
    def get_all_end_effectors(self) -> List[EndEffector]:
        end_effectors = set()
        for phase in self.sequence:
            end_effectors.update(phase.contacts.keys())
        return list(end_effectors)

    # Get the number of end effectors in contact during a phase
    def get_num_stance_during_phase(self, phase_idx: int) -> int:
        count = 0
        for flag in self.sequence[phase_idx].contacts.values():
            count += int(flag)
        return count

    # Get phase from a knot point index
    def get_phase(self, k: int) -> Phase:
        return self.sequence[self.get_phase_idx(k)]

    # Return whether an end effector is in contact at a certain knot point
    def is_in_contact(self, end_effector: EndEffector, k: int) -> bool:
        phase = self.get_phase(k)
        return phase.contacts.get(end_effector, False)

    # Get the frame id from an end effector object (used with pinocchio)
    def get_frame_id(self, end_effector: EndEffector) -> int:
        return end_effector.frame_id

    # Get the flag which specifies whether the end effector has 6 DOF (True) or 3 DOF
    def get_contact_type(self, end_effector: EndEffector) -> bool:
        return end_effector.type_6D

    # Convert the contact type flag into a number (True = 6, False = 3)
    # Used when determining input variable size with casadi
    def get_contact_size(self, end_effector: EndEffector) -> int:
        return 3 * int(end_effector.type_6D) + 3

    # Get the total number of phases in the sequence
    def get_num_phases(self) -> int:
        return len(self.sequence)

    # Convert a knot point into a phase index
    def get_phase_idx(self, k: int) -> int:
        for idx, cum_knot in enumerate(self.cumulative_knots):
            if k < cum_knot:
                return idx
        return None

    def get_time_vec(self, dts: defaultdict(list)):
        t = defaultdict(list)
        k = 0
        tc = 0.0
        for i in range(self.get_num_phases()):
            phase = self.sequence[i]
            t[phase.phase_name] = []
            for ki in range(phase.knot_points):
                t[phase.phase_name].append(tc + dts[i])
                tc += dts[i]
                k += 1
        return t


# Class for setting up pseudospectral collocation
class PseudoSpectralCollocation:
    def __init__(s, degree: int):
        # Degree of interpolating polynomial
        s.degree = degree

        # Get collocation points
        s.tau_root = np.append(0, ca.collocation_points(s.degree, "radau"))

        # Coefficients of the collocation equation
        s.C = np.zeros((s.degree + 1, s.degree + 1))

        # Coefficients of the continuity equation
        s.D = np.zeros(s.degree + 1)

        # Coefficients of the quadrature function
        s.B = np.zeros(s.degree + 1)

        # Construct polynomial basis
        for j in range(s.degree + 1):
            # Construct Lagrange polynomials to get the polynomial basis at the collocation point
            p = np.poly1d([1])
            for r in range(s.degree + 1):
                if r != j:
                    p *= np.poly1d([1, -s.tau_root[r]]) / (s.tau_root[j] - s.tau_root[r])

            # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
            s.D[j] = p(1.0)

            # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the
            # continuity equation
            p_der = np.polyder(p)
            for r in range(s.degree + 1):
                s.C[j, r] = p_der(s.tau_root[r])

            # Evaluate the integral of the polynomial to get the coefficients of the quadrature function
            pint = np.polyint(p)
            s.B[j] = pint(1.0)


# Simple parameters storage
class Parameters:
    def __init__(s, q0: np.ndarray):
        # The initial configration
        s.q0 = q0
        # Gravity vector
        s.g = np.array([0, 0, 9.81])
        # Friction coefficient
        s.mu = 0.7
        # Max swing height
        s.swing_height = 0.15
        # Frequency to interpolate the solution to
        s.desired_frequency = 0.01
        # Problem options for casadi
        s.p_opts = {"expand": True}
        # Solver options for casadi
        s.s_opts = {"max_iter": 200, 
                    # "linear_solver": "ma97"
                    }


# Main optimization class
class CentroidalTrajOpt:
    def __init__(
        s,
        model,
        data,
        viz: MeshcatVisualizer,
        params: Parameters,
        contact_sequence: ContactSequence,
        collocation: PseudoSpectralCollocation,
    ):
        s.model = model
        s.data = data
        s.viz = viz
        s.p = params
        s.cons = contact_sequence
        s.col = collocation
        s.q0 = s.p.q0
        pin.computeTotalMass(s.model, s.data)
        pin.framesForwardKinematics(s.model, s.data, s.q0)

        s.Mtarget = pin.SE3(
            s.data.oMf[s.model.getFrameId("r_foot_v_ft_link")].rotation.copy(),
            np.array([0.0775, 0.2, 0.0]),
        )

        s.viz.addBox("world/box", [0.05, 0.1, 0.2], [1.0, 0.2, 0.2, 0.5])
        s.viz.addBox("world/blue", [0.08] * 3, [0.2, 0.2, 1.0, 0.5])

        for ee in s.cons.get_all_end_effectors():
            s.viz.addSphere(ee.frame_name, [0.07], [0.8, 0.8, 0.2, 0.5])

        s.display_scene(s.q0)

        s.m = s.data.mass[0]
        s.N = s.cons.cumulative_knots[-1]

        s.cmodel = cpin.Model(s.model)
        s.cdata = s.cmodel.createData()

        # Size of joint position vector
        s.nq = s.model.nq
        # Size of joint tangent vector
        s.nv = s.model.nv
        # 3 linear, 3 angular
        s.nm = 3
        # Size of centroidal momentum vector
        s.nh = 2 * s.nm
        # Size of centroidal momentum time derivative vector
        s.ndh = 2 * s.nh

        # Size of the state vector
        s.nx = 2 * s.nh + s.nq + s.nv
        # Size of the state delta vector
        s.ndx = 2 * s.nh + 2 * s.nv

        # State vector where x = [h; hdot; q; v]
        s.cx = ca.SX.sym("x", s.nx, 1)
        # State delta vector
        s.cdx = ca.SX.sym("dx", s.ndx, 1)

        # Contact force and contact wrench
        s.cu = ca.SX.sym("u", 6, 1)
        # Joint tangents are one of the decision variables
        s.cvju = ca.SX.sym("vju", s.nv - 6, 1)

        # Momenta: nh x 1
        s.ch = s.cx[: s.nh]
        # Momenta delta: nh x 1
        s.ch_d = s.cdx[: s.nh]

        # Momenta time derivative: nh x 1
        s.cdh = s.cx[s.nh : s.ndh]
        # Momentum time derivative delta: nh x 1
        s.cdh_d = s.cdx[s.nh : s.ndh]

        # q: nq x 1
        s.cq = s.cx[s.ndh : s.nq + s.ndh]
        # q delta: nv x 1
        s.cq_d = s.cdx[s.ndh : s.nv + s.ndh]

        # qj: (nq - 7) x 1
        s.cqj = s.cq[7:]

        # v: nv x 1
        s.cv = s.cx[s.ndh + s.nq :]
        # v delta: nv x 1
        s.cv_d = s.cdx[s.ndh + s.nv :]

        # v_j: (nv - 6) x 1
        s.cvj = s.cv[6:]

        # f: 3 x 1
        s.cf = s.cu[:3]
        # tau: 3 x 1
        s.ctau = s.cu[3:]

        # Polynomial variables
        s.cxp = ca.SX.sym("xp", s.nx, 1)
        # Polynomial momenta: nh x 1
        s.chp = s.cxp[: s.nh]
        # Polynomial momentum time derivative: nh x 1
        s.cdhp = s.cxp[s.nh : s.ndh]
        # Polynomial q: nq x 1
        s.cqp = s.cxp[s.ndh : s.ndh + s.nq]
        # Polynomial v: nv x 1
        s.cvp = s.cxp[s.ndh + s.nq :]

        # Time helper variable
        s.cdt = ca.SX.sym("dt", 1, 1)

    def compute_casadi_graphs(s):
        # Set up the sym graph to get center of mass quantities
        cpin.centerOfMass(s.cmodel, s.cdata, s.cq, False)
        # Compute the centroidal map outlined in "Centroidal dynamics of a humanoid robot"
        cpin.computeCentroidalMap(s.cmodel, s.cdata, s.cq)
        # Forward kinematics sym graph
        cpin.forwardKinematics(s.cmodel, s.cdata, s.cq, s.cv)
        cpin.updateFramePlacements(s.cmodel, s.cdata)

        # Integrate over the state and the Lie Group with a small delta
        s.cintegrate = ca.Function(
            "integrate",
            [s.cx, s.cdx, s.cdt],
            [
                ca.vertcat(
                    s.ch_d,
                    s.cdh_d,
                    cpin.integrate(s.cmodel, s.cq, s.cq_d * s.cdt),
                    s.cv_d,
                )
            ],
            ["x", "dx", "dt"],
            ["integral(x, dx, dt)"],
        )

        # Centroidal dynamics in continuous time (recall that x = [h; hdot; q; v])
        Ag = s.cdata.Ag
        s.xdot = ca.Function(
            "centroidal_dynamics",
            [s.cx, s.cu, s.cvju],
            [
                ca.vertcat(
                    s.cdh,
                    (s.cf - s.m * s.p.g) / s.m,
                    (s.ctau) / s.m,
                    s.cv,
                    ca.mtimes(
                        ca.inv(Ag[:, :6]),
                        (s.m * s.ch - ca.mtimes(Ag[:, 6:], s.cvju)),
                    ),
                    s.cvju,
                )
            ],
            ["x", "u", "vju"],
            ["xdot"],
        )

        # Running cost
        s.L = ca.Function(
            "running_cost",
            [s.cx, s.cu, s.cvju],
            [
                1e-3 * ca.sumsqr(s.cvju)
                + 1e-4 * ca.sumsqr(s.cf)
                + 1e-4 * ca.sumsqr(s.ctau)
                + 1e1 * ca.sumsqr(s.cqj - s.q0[7:])
            ],
            ["x", "u", "vju"],
            ["L"],
        )
        # # print ("cq: ", s.cq)
        # # print("cv: ", s.cv)
        # # print("cq_d: ", s.cq_d)
        # # print("integrate: ", cpin.integrate(s.cmodel, s.cq, s.cq_d * s.cdt))
        # # print("F: ", s.xdot.sx_out(0))
        # q1 = np.zeros(s.nq)
        # q1[:7] = np.array([0,0,1.0627,0,0,0,1.])
        # # q1 = pin.randomConfiguration(s.model, -np.ones(s.nq), np.ones(s.nq))
        # # q1 = pin.normalize(s.model, q1)
        # # q1[7:] = np.zeros(s.nq - 7)
        # dq1 = np.zeros(s.nv)
        # # dq1 = 2. * np.random.rand(s.nv) - 1.
        # dq1[:6] = np.array([0,0,0,0,0,0])
        # # dq1[6:] = np.zeros(s.nv - 6)

        # print("q0test[0] =  ", q1[0], "; q0test[1] = ", q1[1], "; q0test[2] = ", q1[2], "; q0test[3] = ", q1[3], "; q0test[4] = ", q1[4], "; q0test[5] = ", q1[5], "; q0test[6] = ", q1[6], ";")
        # print("dqtest[0] =  ", dq1[0], "; dqtest[1] = ", dq1[1], "; dqtest[2] = ", dq1[2], "; dqtest[3] = ", dq1[3], "; dqtest[4] = ", dq1[4], "; dqtest[5] = ", dq1[5], ";")
        # inttmp = pin.integrate(s.model, q1, dq1)
        # print("pq0res[0] = ", inttmp[0], "; pq0res[1] = ", inttmp[1], "; pq0res[2] = ", inttmp[2], "; pq0res[3] = ", inttmp[3], "; pq0res[4] = ", inttmp[4], "; pq0res[5] = ", inttmp[5], "; pq0res[6] = ", inttmp[6], ";")

        # Error in contact position
        s.dpcontacts = {}
        # Error in contact velocity
        s.vcontacts = {}

        for ee in s.cons.get_all_end_effectors():
            if not ee.type_6D:
                p0 = s.data.oMf[ee.frame_id].translation.copy()
                s.dpcontacts[ee.frame_name] = ca.Function(
                    f"dpcontact_{ee.frame_name}",
                    [s.cx],
                    [-(s.cdata.oMf[ee.frame_id].inverse().act(ca.SX(p0)))],
                    ["x"],
                    ["3D_pos_error"],
                )
                s.vcontacts[ee.frame_name] = ca.Function(
                    f"vcontact_{ee.frame_name}",
                    [s.cx],
                    [
                        cpin.getFrameVelocity(
                            s.cmodel, s.cdata, ee.frame_id, pin.LOCAL
                        ).linear
                    ],
                    ["x"],
                    ["3D_vel_error"],
                )
            else:
                p0 = s.data.oMf[ee.frame_id]
                s.dpcontacts[ee.frame_name] = ca.Function(
                    f"dpcontact_{ee.frame_name}",
                    [s.cx],
                    [np.zeros(6)],
                    ["x"],
                    ["6D_pos_error"],
                )
                s.vcontacts[ee.frame_name] = ca.Function(
                    f"vcontact_{ee.frame_name}",
                    [s.cx],
                    [
                        cpin.getFrameVelocity(
                            s.cmodel, s.cdata, ee.frame_id, pin.LOCAL
                        ).vector
                    ],
                    [""],
                    ["6D_vel_error"],
                )

        # Get initial contact position
        pin.framesForwardKinematics(s.model, s.data, s.q0)

        # Helper function to get past casadi MX->SX syntax
        s.tau_cross = {
            ee.frame_name: ca.Function(
                f"f_tau_{ee.frame_name}",
                [s.cx, s.cu],
                [
                    ca.cross(
                        s.cdata.oMf[ee.frame_id].translation - s.cdata.com[0],
                        s.cf,
                    )
                ],
                ["x", "u"],
                ["cross(r_c, f)"],
            )
            for ee in s.cons.get_all_end_effectors()
        }

        # Final configuration cost (should only be active for one end effector)
        s.M = {
            ee.frame_name: ca.Function(
                f"terminal_cost_{ee.frame_name}",
                [s.cx],
                [
                    1e4
                    * ca.sumsqr(
                        cpin.log6(
                            s.cdata.oMf[ee.frame_id].inverse() * cpin.SE3(s.Mtarget)
                        ).vector
                    )
                ],
                ["x"],
                ["M"],
            )
            for ee in s.cons.get_all_end_effectors()
        }

        s.z_foot = {
            ee.frame_name: ca.Function(
                f"z_foot_{ee.frame_name}",
                [s.cx],
                [s.cdata.oMf[ee.frame_id].translation[2]],
                ["x"],
                ["z_foot"],
            )
            for ee in s.cons.get_all_end_effectors()
        }

    # Setup the optimization problem using casadi
    def setup_problem(s):
        opti = ca.Opti()
        x0 = np.concatenate([np.zeros(s.ndh), s.q0, np.zeros(s.nv)])
        # Small deviations from the initial state
        # var_dxs is an s.N + 1 size list with s.ndx sized vectors
        var_dxs = [opti.variable(s.ndx) for k in range(s.N + 1)]
        # var_dts holds the contact timing optimization variables
        # (1 per phase, representing the phase duration)
        var_dts = []
        # dts maps the phase to either a float or casadi variable which represents the phase duration
        dts = {}
        # var_us holds the input optimization variables
        # (each element is of size s.cons.get_contact_size(ee))
        var_us = []
        # us maps an end effector frame name to a list which holds all values of u for each knot point
        # This list can contain either casadi variables or numpy zeros,
        # depending on if the end effector is in contact at a certain knot point
        us = defaultdict(list)
        # The centroidal dynamics model assumes we have 'sufficient control authority' and thus we do not care about the joint torques-
        # we consider the joint velocities (var_js) to be the inputs to the system, since with 'sufficient control authority' we can always map any
        # joint acceleration to an associated joint torque using inverse dynamics
        var_vjs = [opti.variable(s.nv - 6) for k in range(s.N)]
        # var_xs is the state integrated from an initial state of np.concatenate([np.zeros(s.ndh), s.q0, np.zeros(s.nv)]) by var_dxs
        var_xs = [s.cintegrate(x0, var_dx, 1) for var_dx in var_dxs]
        # Constrain initial stait
        opti.subject_to(var_xs[0] == x0)

        # Set up cost function
        totalcost = 0
        # Iterate over end effectors to get the appropriate us (numpy zeros or casadi decision variables)
        # depending on if the end effector is in contact at each knot point
        for ee in s.cons.get_all_end_effectors():
            # Initialize the list
            us[ee.frame_name] = []
            for k in range(s.N):
                # 6 DOF contact or 3 DOF contact
                legsize = s.cons.get_contact_size(ee)
                if s.cons.is_in_contact(ee, k):
                    u = opti.variable(legsize)
                    var_us.append(u)
                    us[ee.frame_name].append(u)
                else:
                    us[ee.frame_name].append(np.zeros(legsize))

        # Iterate over the phases to get the appropriate dt (a float or casadi decision variable)
        # depending on if the phase.timing_var flag is True
        for phase in s.cons.sequence:
            if phase.timing_var:
                dt_var = opti.variable(1)
                var_dts.append(dt_var)
                dts[phase.phase_name] = dt_var
            else:
                dts[phase.phase_name] = phase.fixed_timing

        # "Lift" initial conditions
        x_k = var_xs[0]
        dx_k = var_dxs[0]
        # Iterate over each knot point
        for k in range(s.N):
            # Get the relevant phase and dt
            phase = s.cons.get_phase(k)
            i = s.cons.get_phase_idx(k)
            dt = dts[phase.phase_name]

            # grf is the total ground reaction force acting in the COM frame aligned with the inertial frame
            grf = np.zeros(3)
            # tau is the total contact wrench in the COM frame aligned with the inertial frame
            tau = np.zeros(3)
            # Iterate over all end effectors
            for ee in s.cons.get_all_end_effectors():
                # Get the stacked force and torque of the end effector at the current knot point
                # F is a vector of zeros if the end effector is not in contact,
                # and is a casadi decision variable if it is in contact
                F = us[ee.frame_name][k]
                # Check if the end effector is in contact at the current knot point
                if s.cons.is_in_contact(ee, k):
                    # Add the contribution of the forces and torques in the COM frame aligned with the inertial frame
                    grf += F[:3]
                    tau += s.tau_cross[ee.frame_name](x_k, F)
                    # Add all the contact constraints:
                    # Zero velocity at the contact point
                    opti.subject_to(s.vcontacts[ee.frame_name](x_k) == 0)
                    # Friction cone approximation constraint
                    opti.subject_to(F[0] / F[2] > -s.p.mu)
                    opti.subject_to(F[0] / F[2] < s.p.mu)
                    opti.subject_to(F[1] / F[2] > -s.p.mu)
                    opti.subject_to(F[1] / F[2] < s.p.mu)
                    # Unilateral contact constraint
                    # (normal force must be positive- we can push against the ground but we can't pull!)
                    opti.subject_to(F[2] >= 0)
                    # Set the initial guess for the force variable
                    # (an initial guess of 0 is infeasible because of the friction cone constraints)
                    opti.set_initial(
                        F[2], s.m * s.p.g[2] / s.cons.get_num_stance_during_phase(i)
                    )
                    # I think this constraint is wrong! F[3:] is the wrench at the point of contact,
                    # but tau is te total contact wrench in the COM frame aligned with the inertial frame
                    if s.cons.get_contact_size(ee) == 6:
                        tau += F[3:]
                else:
                    kphase = k
                    if i > 0:
                        kphase = k - s.cons.cumulative_knots[i - 1]
                    # Add swing constraints
                    # if kphase == phase.knot_points / 2:
                    #     opti.subject_to(s.swing_p[ee.frame_name](x_k) == 0)
                    # opti.subject_to(s.swing_v[ee.frame_name](x_k) == 0)
                    opti.subject_to(s.z_foot[ee.frame_name](x_k) >= 0)
            # Stack the total ground reaction force and contact wrench
            u_k = ca.vertcat(grf, tau)
            # Get the joint velocities at the current knot point
            vj_k = var_vjs[k]

            # State at collocation points (s.col.degree decision variables)
            dx_c = []
            for j in range(s.col.degree):
                dx_kj = opti.variable(s.ndx)
                dx_c.append(dx_kj)

            # Loop over collocation points
            dx_k_end = s.col.D[0] * dx_k
            for j in range(1, s.col.degree + 1):
                dt_j = (s.col.tau_root[j] - s.col.tau_root[j - 1]) * dt
                # Expression for the state derivative at the collocation point
                dx_p = s.col.C[0, j] * dx_k
                for r in range(s.col.degree):
                    dx_p += s.col.C[r + 1, j] * dx_c[r]

                # Append collocation equations
                x_c = s.cintegrate(x_k, dx_c[j - 1], dt_j)
                fj = s.xdot(x_c, u_k, vj_k)
                qj = s.L(x_c, u_k, vj_k)

                # Since the state includes a quaternion which is in a Lie Group,
                # we need to find the tangent vector which transports us from
                # the quaternion at x_k to -x_p, where x_p is the approximated
                # state derivative at the collocation point. However, I'm not sure that this is right,
                # since the 'expression for the state derivative at the collocation point' is in terms
                # of the 4 x 1 vector defining a quaternion, but the state derivative of the quaternion
                # is actually a 3 x 1 tangent vector. How do you rectify this?
                # dx_col = s.dx_col(x_k, x_p)

                # By multiplying fj with h, we're effectively rescaling the state's rate of change from the
                # time scale of the problem to the normalized time scale used by Lagrange polynomials
                opti.subject_to(dt * fj - dx_p == 0)
                # Add contribution to the end state
                dx_k_end += s.col.D[j] * dx_c[j - 1]
                # Add contribution to quadrature function
                totalcost += s.col.B[j] * qj * dt

            # Get the NLP variable for the state at end of interval
            x_k = var_xs[k + 1]
            dx_k = var_dxs[k + 1]
            # Add equality constraint
            opti.subject_to(dx_k_end - dx_k == 0)

        totalcost += s.M["r_foot_v_ft_link"](x_k)

        get_dts = ca.Function(
            "get_dts",
            [ca.vertcat(*var_dts)],
            [ca.vertcat(*[dts[phase.phase_name] for phase in s.cons.sequence])],
        )

        print("Optimizing...")

        # SOLVE
        opti.minimize(totalcost)
        # Set numerical backends
        opti.solver("ipopt", s.p.p_opts, s.p.s_opts)
        opti.callback(
            lambda i: s.display_scene(
                opti.debug.value(var_xs[-1][s.ndh : s.ndh + s.nq])
            )
        )

        # Caution: in case the solver does not converge, we are picking the candidate values
        # at the last iteration in opti.debug, and they are NO guarantee of what they mean.
        try:
            sol = opti.solve_limited()
            sol_xs = [opti.value(var_x) for var_x in var_xs]
            sol_dts = get_dts(ca.vertcat(*[opti.value(var_dt) for var_dt in var_dts]))
        except:
            print("ERROR in convergence, plotting debug info.")
            sol_xs = [opti.debug.value(var_x) for var_x in var_xs]
            sol_dts = get_dts(
                ca.vertcat(*[opti.debug.value(var_dt) for var_dt in var_dts])
            )
        print("***** Display the resulting trajectory ...")
        xdes = np.array([x[s.ndh : s.ndh + s.nq] for x in sol_xs])
        sol_dts = sol_dts.full().flatten()

        t = s.cons.get_time_vec(sol_dts)
        xnew, tnew = s.interpolate(xdes, t)

        while True:
            s.display_traj(xnew)
            xnewrev = np.fliplr(xnew)
            s.display_traj(xnewrev)

    def interpolate(s, x, t) -> (np.ndarray, np.ndarray):
        cumknots = 0
        xnew = np.array([])
        tnew = np.array([])
        flag = True
        for phase in s.cons.sequence:
            tphase = t[phase.phase_name]
            xphase = x[cumknots : cumknots + phase.knot_points]
            tgrid_new = np.linspace(
                tphase[0],
                tphase[-1],
                int((tphase[-1] - tphase[0]) / s.p.desired_frequency),
            )
            x_phase_new = np.array([])
            for i in range(s.nq):
                xi_phase_opt_new = np.array([])
                if i < 3 or i > 6:
                    f = interp1d(tphase, xphase[:, i], kind="cubic")
                    xi_phase_opt_new = f(tgrid_new)
                if i == 3:
                    slerp = Slerp(tphase, R.from_quat(xphase[:, 3:7]))
                    xi_phase_opt_new = np.transpose(slerp(tgrid_new).as_quat())
                if i == 0:
                    x_phase_new = xi_phase_opt_new
                elif i <= 3 or i > 6:
                    x_phase_new = np.vstack((x_phase_new, xi_phase_opt_new))
            cumknots += phase.knot_points
            if flag:
                xnew = x_phase_new
                tnew = tgrid_new
                flag = False
            else:
                xnew = np.hstack((xnew, x_phase_new))
                tnew = np.hstack((tnew, tgrid_new))
        return xnew, tnew

    def display_scene(s, q: np.ndarray, dt=1e-1):
        """
        Given the robot configuration, display:
        - the robot
        - a box representing endEffector_ID
        - a box representing Mtarget
        """
        pin.framesForwardKinematics(s.model, s.data, q)
        M = s.data.oMf[s.model.getFrameId("r_foot_v_ft_link")]
        s.viz.applyConfiguration("world/box", s.Mtarget)
        s.viz.applyConfiguration("world/blue", M)
        for ee in s.cons.get_all_end_effectors():
            s.viz.applyConfiguration(ee.frame_name, s.data.oMf[ee.frame_id])
        s.viz.display(q)
        time.sleep(dt)

    def display_traj(s, qs: np.ndarray):
        for k in range(np.size(qs, 1)):
            s.display_scene(qs[:, k], s.p.desired_frequency)


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
        1,
        0,
        0,
        0,
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

viz = MeshcatVisualizer(robot)
viz.display(robot.q0)

params = Parameters(robot.q0)

l_name = "l_foot_v_ft_link"
r_name = "r_foot_v_ft_link"

# Define two end effectors
ee_left_foot = EndEffector(
    frame_name=l_name, frame_id=model.getFrameId(l_name), type_6D=True
)
ee_right_foot = EndEffector(
    frame_name=r_name, frame_id=model.getFrameId(r_name), type_6D=True
)
phase_1_period = 0.25
phase_1_knot_points = 25
phase_1_fixed_timing = phase_1_period / phase_1_knot_points
contacts_phase1 = frozendict({ee_left_foot: True, ee_right_foot: True})
phase1 = Phase(
    contacts=contacts_phase1,
    phase_name="phase_1",
    fixed_timing=phase_1_fixed_timing,
    knot_points=phase_1_knot_points,
)
phase_2_period = 0.25
phase_2_knot_points = 25
phase_2_fixed_timing = phase_2_period / phase_2_knot_points
contacts_phase2 = frozendict({ee_left_foot: True, ee_right_foot: False})
phase2 = Phase(
    contacts=contacts_phase2,
    phase_name="phase_2",
    fixed_timing=phase_2_fixed_timing,
    knot_points=phase_2_knot_points,
)
phase_3_period = 0.01
phase_3_knot_points = 5
phase_3_fixed_timing = phase_3_period / phase_3_knot_points
contacts_phase3 = frozendict({ee_left_foot: True, ee_right_foot: True})
phase3 = Phase(
    contacts=contacts_phase3,
    phase_name="phase_3",
    fixed_timing=phase_3_fixed_timing,
    knot_points=phase_3_knot_points,
)

# Create a contact sequence and add phases
contact_seq = ContactSequence()
contact_seq.add_phase(phase1)
contact_seq.add_phase(phase2)
# contact_seq.add_phase(phase3)


collocation = PseudoSpectralCollocation(1)

traj_opt = CentroidalTrajOpt(model, data, viz, params, contact_seq, collocation)
traj_opt.compute_casadi_graphs()
print("Finished computing casadi graphs")
traj_opt.setup_problem()
