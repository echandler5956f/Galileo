#!/home/quant/ros_ws/src/huron_centroidal/.venv/bin/python3

import rospy
from std_msgs.msg import Float64MultiArray
from sensor_msgs.msg import JointState
from std_srvs.srv import Empty, EmptyRequest
from gazebo_msgs.srv import GetModelState, GetModelStateRequest

from types import SimpleNamespace


import time
import casadi
import numpy as np
import pinocchio as pin
from pinocchio import casadi as cpin
from pinocchio.robot_wrapper import RobotWrapper
import numpy as np
import casadi as ca

import rospy
from gazebo_msgs.srv import SpawnModel


xc = np.zeros((19, 1))
vc = np.zeros((18, 1))


def joint_callback(data):
# 0 l_ankle_pitch_joint
# 1 l_ankle_roll_joint
# 2 l_hip_pitch_joint
# 3 l_hip_roll_joint
# 4 l_hip_yaw_joint
# 5 l_knee_pitch_joint
# 6 r_ankle_pitch_joint
# 7 r_ankle_roll_joint
# 8 r_hip_pitch_joint
# 9 r_hip_roll_joint
# 10 r_hip_yaw_joint
# 11 r_knee_pitch_joint
    global xc, vc, tauc
    pos = data.position
    vel = data.velocity
    xc[0] = pos[4]
    xc[1] = pos[3]
    xc[2] = pos[2]
    xc[3] = pos[5]
    xc[4] = pos[0]
    xc[5] = pos[1]

    xc[6] = pos[10]
    xc[7] = pos[9]
    xc[8] = pos[8]
    xc[9] = pos[11]
    xc[10] = pos[6]
    xc[11] = pos[7]

    vc[0] = vel[4]
    vc[1] = vel[3]
    vc[2] = vel[2]
    vc[3] = vel[5]
    vc[4] = vel[0]
    vc[5] = vel[1]

    vc[6] = vel[10]
    vc[7] = vel[9]
    vc[8] = vel[8]
    vc[9] = pos[11]
    vc[10] = vel[6]
    vc[11] = vel[7]


class SimpleInterface:
    def __init__(self) -> None:
        # --- Load robot model
        builder = RobotWrapper.BuildFromURDF
        self.robot = builder(
            "/home/quant/ros_ws/src/HURON-Model/huron_description/urdf/huron_cheat.urdf",
            ["/home/quant/ros_ws/src/HURON-Model/huron_description/urdf"],
            None,
        )
        self.robot.q0 = np.array(
        [
            0,
            0,
            1.1064,
            0,
            0,
            0,
            1,
            0.0000,
            0.0000,
            -0.1672,
            0.3922,
            -0.2251,
            0.0000,
            0.0000,
            0.0000,
            -0.1672,
            0.3922,
            -0.2251,
            0.0000,
        ]
        )

        self.model = self.robot.model
        self.data = self.model.createData()
        pin.framesForwardKinematics(self.model, self.data, self.robot.q0)
        pin.computeTotalMass(self.model, self.data)
        pin.updateFramePlacements(self.model, self.data)
        self.mass = self.data.mass[0]

    def get_com(self, q):
        pin.centerOfMass(self.model, self.data, q, False)
        return self.data.com[0]

    def get_vcom(self, q, v):
        pin.centerOfMass(self.model, self.data, q, v, False)
        return self.data.vcom[0]

    def get_A(self, q):
        pin.computeCentroidalMap(self.model, self.data, q)
        return self.data.Ag

    def get_Adot(self, q, v):
        pin.computeCentroidalMapTimeVariation(self.model, self.data, q, v)
        return self.data.dAg

    def get_h(self, q, v):
        return np.matmul(self.get_A(q), v)

    def get_hdot(self, q, v, a):
        return np.matmul(self.get_A(q), a) + np.matmul(self.get_Adot(q, v), v)

    def get_left_foot_pos(self, q):
        pin.forwardKinematics(self.model, self.data, q)
        pin.updateFramePlacements(self.model, self.data)
        return self.data.oMf[self.model.getFrameId("l_foot_v_ft_link")].translation

    def get_right_foot_pos(self, q):
        pin.forwardKinematics(self.model, self.data, q)
        pin.updateFramePlacements(self.model, self.data)
        return self.data.oMf[self.model.getFrameId("r_foot_v_ft_link")].translation

    def compute_dynamic_components(self, q, qdot):
        # Mass matrix
        pin.crba(self.model, self.data, q)
        pin.nonLinearEffects(self.model, self.data, q, qdot)
        M_symbolic = self.data.M
        # mass_matrix_func = ca.Function("mass_matrix_func", [q], [M_symbolic])

        # Coriolis matrix
        C_symbolic = ca.mtimes(self.data.C, qdot)
        # coriolis_matrix_func = ca.Function("coriolis_matrix_func", [q, qdot], [C_symbolic])

        # Gravity vector
        g_sym = self.data.g
        # gravity_vector_func = ca.Function("gravity_vector", [q], [g_sym])

        # Torque selection matrix
        S = np.zeros((self.model.nv, self.model.nv))
        # 12 x 18
        S[6:, 6:] = np.eye(self.model.nv - 6)

        return M_symbolic, C_symbolic, g_sym, S


def main():
    global xc, vc, tauc

    q0 = np.reshape(np.array(
        [
            0,
            0,
            1.1064,
            0,
            0,
            0,
            1,
            0.0000,
            0.0000,
            -0.1672,
            0.3922,
            -0.2251,
            0.0000,
            0.0000,
            0.0000,
            -0.1672,
            0.3922,
            -0.2251,
            0.0000,
        ]
    ), (19, 1))

    xc = q0
    vc = np.zeros((18, 1))

    print("Start SimpleInterface setup")

    ibrahim = SimpleInterface()

    model = ibrahim.model
    data = ibrahim.data
    mass = ibrahim.mass

    print("SimpleInterface setup done")

    # rospy.wait_for_service('/gazebo/unpause_physics')
    # rospy.wait_for_service('/gazebo/pause_physics')
    rospy.wait_for_service("/gazebo/get_model_state")
    unpause_physics_client = rospy.ServiceProxy("/gazebo/unpause_physics", Empty)
    # pause_physics_client = rospy.ServiceProxy('/gazebo/pause_physics', Empty)

    get_model_state_client = rospy.ServiceProxy(
        "/gazebo/get_model_state", GetModelState
    )
    rospy.wait_for_service("/gazebo/spawn_urdf_model")
    spawn_model_service = rospy.ServiceProxy("/gazebo/spawn_urdf_model", SpawnModel)

    rospy.Subscriber("/huron/joint_states", JointState, joint_callback, queue_size=100)
    tau_publisher = rospy.Publisher(
        "/huron/joint_group_effort_controller/command", Float64MultiArray, queue_size=100
    )
    rospy.init_node("ibrahim_test", anonymous=True)

    rospy.wait_for_service("/gazebo/unpause_physics")

    # or write the standing control as follows.......................................
    l_name = "l_foot_v_ft_link"
    r_name = "r_foot_v_ft_link"

    l_ID = model.getFrameId(l_name)
    r_ID = model.getFrameId(r_name)

    # initialize symbolic model and data
    cmodel = cpin.Model(model)
    cdata = cmodel.createData()
    # number of position variables (19)
    nq = model.nq
    # number of velocity variables (18)
    nv = model.nv

    # build the expresion graphs (note- cq and cdq are not decision variables
    # they are merely used to build the symbolic functions which you can plug numeric values into)
    cq0 = ca.SX.sym("q0", nq, 1)
    cq = ca.SX.sym("q", nq, 1)
    cdq = ca.SX.sym("dq", nv, 1)
    cddq = ca.SX.sym("ddq", nv, 1)
    cdt = ca.SX.sym("dt", 1, 1)

    # Define symbolic variables for left and right contact forces
    f_left = ca.SX.sym("f_left", 6, 1)  # Assuming a 3D force for left contact
    f_right = ca.SX.sym("f_right", 6, 1)  # Assuming a 3D force for right contact

    # Define symbolic tau
    tau = ca.SX.sym("tau", nv, 1)

    # blah blah
    cpin.centerOfMass(cmodel, cdata, cq, cdq, False)
    cpin.computeCentroidalMap(cmodel, cdata, cq)
    cpin.computeCentroidalMapTimeVariation(cmodel, cdata, cq, cdq)
    cpin.forwardKinematics(cmodel, cdata, cq, cdq, cddq)
    cpin.updateFramePlacements(cmodel, cdata)

    # get CMM and time derivative
    Ag = cdata.Ag
    dAg = cdata.dAg

    # get contact jacobians and derivatives in the proper frame
    J1 = cpin.computeFrameJacobian(cmodel, cdata, cq, l_ID, pin.LOCAL)
    J2 = cpin.computeFrameJacobian(cmodel, cdata, cq, r_ID, pin.LOCAL)
    J1dot = cpin.frameJacobianTimeVariation(cmodel, cdata, cq, cdq, l_ID, pin.LOCAL)
    J2dot = cpin.frameJacobianTimeVariation(cmodel, cdata, cq, cdq, r_ID, pin.LOCAL)

    # # get the spatial acceleration of the contacts
    acc1 = cpin.getFrameAcceleration(cmodel, cdata, l_ID, pin.LOCAL).vector
    acc2 = cpin.getFrameAcceleration(cmodel, cdata, r_ID, pin.LOCAL).vector

    vel1 = cpin.getFrameVelocity(cmodel, cdata, l_ID, pin.LOCAL).vector
    vel2 = cpin.getFrameVelocity(cmodel, cdata, r_ID, pin.LOCAL).vector
    print(vel1.shape)
    print(J1.shape)

    # tuning weights
    w1 = np.eye(3)
    w2 = np.eye(3)
    w3 = np.eye(nv)
    beta1 = 1.0
    beta2 = 1.0
    beta3 = 1.0

    # function for the symbolic angular momentum rate of change
    kdot_fun = ca.Function(
        "ang_mom_dot",
        [cq, cdq, cddq],
        [ca.mtimes(Ag[:3, :], cddq) + ca.mtimes(dAg[:3, :], cdq)],
    )
    # function for the symbolic linear momentum rate of change
    ldot_fun = ca.Function(
        "lin_mom_dot",
        [cq, cdq, cddq],
        [ca.mtimes(Ag[3:, :], cddq) + ca.mtimes(dAg[3:, :], cdq)],
    )

    # enforces the constraint in (Eq. 13) for both contacts (note that this ONLY enforces the relationship between contact acceleration and joint acceleration. nothing to guarantee unilateral forces)
    acc_constraint = ca.Function(
        "acc_constraint",
        [cq, cdq, cddq],
        [
            ca.vertcat(
                # ca.mtimes(J1, cddq) + ca.mtimes(J1dot, cdq),# - ca.mtimes(0. * np.zeros(6), ca.mtimes(J1, cdq)),
                # ca.mtimes(J2, cddq) + ca.mtimes(J2dot, cdq) #- ca.mtimes(0. * np.zeros(6), ca.mtimes(J2, cdq))
                acc1,
                acc2,
            )
        ],
    )
    DT = 1e-3
    cintegrate = ca.Function(
        "integrate",
        [cq, cdq, cddq, cdt],
        [cpin.integrate(cmodel, cq, cdq * cdt + cddq * cdt * cdt)],
    )

    cdifference = ca.Function(
        "difference",
        [cq0, cq, cdq, cddq, cdt],
        [cpin.difference(cmodel, cq0, cintegrate(cq, cdq, cddq, cdt))],
    )
    vcontact_l = ca.Function(
        "vcontact_l",
        [cq, cdq],
        [cpin.getFrameVelocity(cmodel, cdata, l_ID, pin.LOCAL).vector],
    )
    vcontact_r = ca.Function(
        "vcontact_r",
        [cq, cdq],
        [cpin.getFrameVelocity(cmodel, cdata, r_ID, pin.LOCAL).vector],
    )
    xc = q0[7:]
    vc = np.zeros((12, 1))

    i = 0
    unpause_physics_client(EmptyRequest())
    while True:
        seconds = rospy.get_time()
        ms = get_model_state_client(
            GetModelStateRequest("huron_description", "ground_plane")
        )
        p = ms.pose.position
        o = ms.pose.orientation
        vl = ms.twist.linear
        al = ms.twist.angular
        q = np.vstack((np.reshape(np.array([p.x, p.y, p.z, o.x, o.y, o.z, o.w]), (7, 1)), xc))
        v = np.vstack((np.reshape(np.array([vl.x, vl.y, vl.z, al.x, al.y, al.z]), (6, 1)), vc))

        # desired conf
        q_des = q0  # for bending knee configuration
        qdot_des = np.zeros((18, 1))
        # reference acceleration
        qdotdot_ref = np.zeros((18, 1))
        # reference com
        c_r = np.array(
            [[0.0255], [0.0], [0.6244]]
        )  # for bending knee configuration

        # current com
        c = data.com[0]
        # print("com position", c)
        # print("mass=", mass)

        # reference com velocity
        cdot_r = np.zeros((3, 1))
        # current com velocity
        cdot = data.vcom[0]

        # current joint angles
        q_cur = q
        # current joint velocities
        qdot_cur = v

        # tuning variables
        k_l = 1.0
        d_l = 1.0
        k_h = 1.0
        d_h = 1.0
        k_t = 1.0
        d_t = 1.0

        # desired linear momentum rate of change (Eq. 6)
        Ldot_des = k_l * mass * (c_r - c) - d_l * mass * (cdot_r - cdot)

        # desired angular momentum change (Eq. 9)
        Hdot_des = np.zeros((3, 1))
        # desired acceleration (Eq. 11)
        # qdotdot_des = (
        #     k_t * (np.reshape(pin.difference(model, q_des, q_cur), (18, 1)))
        #     + d_t * (qdot_des - qdot_cur)
        #     + qdotdot_ref
        # )
        qdotdot_des = 240.0*np.reshape((q_des[7:] - q_cur[7:]), (12, 1))
        # print(qdotdot_des)

        # print("Start opti")

        # initialize the optimizer
        opti = ca.Opti("conic")
        # create the decision variable (qdotdot)
        var_ddq = opti.variable(nv)
        # create decision variables for tau, f_left, and f_right
        var_tau = opti.variable(nv)
        var_f_left = opti.variable(6)
        var_f_right = opti.variable(6)

        # a scalar objective (Eq. 7, 10, 12)
        objective = 0
        # objective = ca.sumsqr(w3 * (qdotdot_des - var_ddq)) + ca.sumsqr(var_tau)
        objective = objective + (
            beta1 * ca.sumsqr(w1 * (Hdot_des - kdot_fun(q_cur, qdot_cur, var_ddq)))
            + beta2 * ca.sumsqr(w2 * (Ldot_des - ldot_fun(q_cur, qdot_cur, var_ddq)))
            + ca.sumsqr((qdotdot_des - var_ddq[6:]))
            + ca.sumsqr(var_tau)
        )
        # print(objective)

        # activate the acceleration contraint
        opti.subject_to(acc_constraint(q_cur, qdot_cur, var_ddq) == 0)

        J1num = pin.computeFrameJacobian(model, data, q_cur, l_ID, pin.LOCAL)
        J2num = pin.computeFrameJacobian(model, data, q_cur, r_ID, pin.LOCAL)

        # Call compute_dynamic_components function
        M_dyn, C_dyn, g_dyn, S_dyn = ibrahim.compute_dynamic_components(q_cur, qdot_cur)
        g_dyn = ca.reshape(g_dyn, 18, 1)

        A_dyn = ca.vertcat(
            ca.horzcat(
                M_dyn, -S_dyn.transpose(), -J1num.transpose(), -J2num.transpose()
            ),
        )
        b_dyn = -C_dyn - g_dyn  # Right-hand side of the linear equality

        # Add the linear equality constraint
        opti.subject_to(
            ca.mtimes(A_dyn, ca.vertcat(var_ddq, var_tau, var_f_left, var_f_right))
            == b_dyn
        )

        # Add an additional torque constraint (adjust as needed)
        extra_torque_constraint = (
            320  # Replace with your desired torque constraint value
        )
        opti.subject_to(var_tau < extra_torque_constraint)
        opti.subject_to(var_tau > -extra_torque_constraint)

        if i > 0:
            # opti.set_initial(opti._lam_g(), prev_lam_g)
            opti.set_initial(opti.x, prev_x)

        opti.minimize(objective)

        opti.solver("osqp", {"verbose": False})

        sol = opti.solve_limited()
        # prev_lam_g = opti.value(opti._lam_g())
        prev_x = opti.value(opti.x)
        optimized_tau = opti.value(var_tau)
        last_12_torques = optimized_tau[6:]

        marray = Float64MultiArray()
        marray.data = last_12_torques
        tau_publisher.publish(marray)
        print(optimized_tau)
        i = i + 1
        DT = rospy.get_time() - seconds


if __name__ == "__main__":
    try:
        main()
    except rospy.ROSInterruptException:
        pass
