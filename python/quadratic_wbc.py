#!/home/quant/ros_ws/src/huron_centroidal/.venv/bin/python3

import rospy
from std_msgs.msg import Float64MultiArray
from sensor_msgs.msg import JointState
from std_srvs.srv import Empty, EmptyRequest

from gazebo_msgs.srv import (
    GetModelState,
    GetModelStateRequest,
)

import casadi as ca
import numpy as np
import pinocchio as pin
from pinocchio import casadi as cpin
from pinocchio.robot_wrapper import RobotWrapper
import numpy as np

xc = np.zeros((12, 1))
vc = np.zeros((12, 1))


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
    vc[9] = vel[11]
    vc[10] = vel[6]
    vc[11] = vel[7]

def hook():
    global tau_publisher
    marray = Float64MultiArray()
    marray.data = np.zeros((12, 1))
    tau_publisher.publish(marray)

def main():
    global xc, vc, tau_publisher

    q0 = np.reshape(
        np.array(
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
        ),
        (19, 1),
    )

    v0 = np.zeros((18, 1))

    xc = q0[7:]
    rospy.wait_for_service("/gazebo/get_model_state")
    unpause_physics_client = rospy.ServiceProxy("/gazebo/unpause_physics", Empty)
    pause_physics_client = rospy.ServiceProxy("/gazebo/pause_physics", Empty)

    get_model_state_client = rospy.ServiceProxy(
        "/gazebo/get_model_state", GetModelState
    )
    rospy.Subscriber("/huron/joint_states", JointState, joint_callback, queue_size=100)
    tau_publisher = rospy.Publisher(
        "/huron/joint_group_effort_controller/command",
        Float64MultiArray,
        queue_size=100,
    )

    rospy.init_node("ibrahim_test", anonymous=True)
    rospy.wait_for_service("/gazebo/unpause_physics")

    rospy.on_shutdown(hook)

    builder = RobotWrapper.BuildFromURDF
    robot = builder(
        "/home/quant/ros_ws/src/HURON-Model/huron_description/urdf/huron_cheat.urdf",
        ["/home/quant/ros_ws/src/HURON-Model/huron_description/urdf"],
        None,
    )
    robot.q0 = q0
    model = robot.model
    data = model.createData()

    l_name = "l_foot_v_ft_link"
    r_name = "r_foot_v_ft_link"

    l_ID = model.getFrameId(l_name)
    r_ID = model.getFrameId(r_name)

    pin.framesForwardKinematics(model, data, robot.q0)
    pin.computeTotalMass(model, data)
    pin.updateFramePlacements(model, data)

    mass = data.mass[0]

    # number of position variables (19)
    nq = model.nq
    # number of velocity variables (18)
    nv = model.nv

    cmodel = cpin.Model(model)
    cdata = cmodel.createData()

    cq = ca.SX.sym("q", nq, 1)
    cdq = ca.SX.sym("dq", nv, 1)
    cddq = ca.SX.sym("ddq", nv, 1)

    cf_left = ca.SX.sym("f_left", 6, 1)
    cf_right = ca.SX.sym("f_right", 6, 1)
    ctau = ca.SX.sym("tau", nv, 1)

    cpin.centerOfMass(cmodel, cdata, cq, cdq, False)
    cpin.forwardKinematics(cmodel, cdata, cq, cdq, cddq)
    cpin.computeCentroidalMap(cmodel, cdata, cq)
    cpin.computeCentroidalMapTimeVariation(cmodel, cdata, cq, cdq)
    cpin.updateFramePlacements(cmodel, cdata)

    cpin.computeJointJacobians(cmodel, cdata, cq)
    cpin.computeJointJacobiansTimeVariation(cmodel, cdata, cq, cdq)

    # get CMM and time derivative
    Ag = cdata.Ag
    dAg = cdata.dAg

    J = cdata.J
    Jdot = cdata.dJ

    J1 = cpin.computeFrameJacobian(cmodel, cdata, cq, l_ID, pin.LOCAL)
    J2 = cpin.computeFrameJacobian(cmodel, cdata, cq, r_ID, pin.LOCAL)
    J1dot = cpin.frameJacobianTimeVariation(cmodel, cdata, cq, cdq, l_ID, pin.LOCAL)
    J2dot = cpin.frameJacobianTimeVariation(cmodel, cdata, cq, cdq, r_ID, pin.LOCAL)

    acc1 = cpin.getFrameAcceleration(cmodel, cdata, l_ID, pin.LOCAL).vector
    acc2 = cpin.getFrameAcceleration(cmodel, cdata, r_ID, pin.LOCAL).vector

    vel1 = cpin.getFrameVelocity(cmodel, cdata, l_ID, pin.LOCAL).vector
    vel2 = cpin.getFrameVelocity(cmodel, cdata, r_ID, pin.LOCAL).vector

    cpin.crba(cmodel, cdata, cq)
    cpin.computeCoriolisMatrix(cmodel, cdata, cq, cdq)
    cpin.computeGeneralizedGravity(cmodel, cdata, cq)

    Msym = cdata.M
    Csym = ca.mtimes(cdata.C, cdq)
    gsym = cdata.g
    S = np.zeros((nv, model.nv))
    S[6:, 6:] = np.eye(nv - 6)
    # print("M:\n", Msym)
    # print("C:\n", Csym)
    # print("g:\n", gsym)

    A = ca.vertcat(ca.horzcat(Msym, -S.T, -J1.T, J2.T))
    b = -Csym - gsym

    cdynamics = ca.Function(
        "dynamics",
        [
            cq,
            cdq,
            cddq,
            ctau,
            cf_left,
            cf_right,
        ],
        [ca.mtimes(A, ca.vertcat(cddq, ctau, cf_left, cf_right)) - b],
    )

    ccom = ca.Function("com", [cq], [cdata.com[0]])

    cvcom = ca.Function("vcom", [cq, cdq], [cdata.vcom[0]])

    vel1_constraint = ca.Function(
        "vel_constraint",
        [cq, cdq],
        [
            ca.vertcat(
                cpin.getFrameVelocity(cmodel, cdata, l_ID, pin.LOCAL).vector,
            )
        ],
    )

    vel2_constraint = ca.Function(
        "vel_constraint",
        [cq, cdq],
        [
            ca.vertcat(
                cpin.getFrameVelocity(cmodel, cdata, r_ID, pin.LOCAL).vector,
            )
        ],
    )

    acc1_constraint = ca.Function(
        "acc1_constraint",
        [cq, cdq, cddq],
        [
            ca.vertcat(
                cpin.getFrameAcceleration(cmodel, cdata, l_ID, pin.LOCAL).vector,
            )
        ],
    )

    acc2_constraint = ca.Function(
        "acc2_constraint",
        [cq, cdq, cddq],
        [
            ca.vertcat(
                cpin.getFrameAcceleration(cmodel, cdata, r_ID, pin.LOCAL).vector,
            )
        ],
    )

    jac_cost = ca.Function(
        "jac_cost",
        [cq, cdq, cddq],
        [ca.vertcat(ca.mtimes(J, cddq) + ca.mtimes(Jdot, cdq))],
    )

     # function for the symbolic angular momentum rate of change
    Hdot_fun = ca.Function(
        "ang_mom_dot",
        [cq, cdq, cddq],
        [ca.mtimes(Ag[:3, :], cddq) + ca.mtimes(dAg[:3, :], cdq)],
    )
    # function for the symbolic linear momentum rate of change
    Ldot_fun = ca.Function(
        "lin_mom_dot",
        [cq, cdq, cddq],
        [ca.mtimes(Ag[3:, :], cddq) + ca.mtimes(dAg[3:, :], cdq)],
    )

    getM = ca.Function("M", [cq], [ca.vertcat(Msym)])
    getC = ca.Function("C", [cq, cdq], [ca.vertcat(Csym)])
    getG = ca.Function("g", [cq], [ca.vertcat(gsym)])

    getStaticTorque = ca.Function(
        "static_tau",
        [cq, cdq, cf_left, cf_right],
        [ca.vertcat(gsym - ca.mtimes(ca.transpose(J), cf_left + cf_right))],
    )

    # print("M:\n", getM(q0))
    # print("C:\n", getC(q0, v0))
    # print("g:\n", getG(q0))

    mu = 0.7
    num_contacts = 2
    A_i = np.asarray(
        [
            [0, 0, 0, 1, 0, -mu],  # pyramid approximation of CWC for one
            [0, 0, 0, -1, 0, -mu],  # contact force f \in R^3
            [0, 0, 0, 0, 1, -mu],
            [0, 0, 0, 0, -1, -mu],
        ]
    )
    fric_constraint = ca.Function(
        "fric_constraint",
        [cf_left, cf_right],
        [
            ca.vertcat(
                ca.mtimes(
                    np.kron(np.eye(num_contacts), A_i), ca.vertcat(cf_left, cf_right)
                )
            )
        ],
    )
    com_ref = ccom(q0)
    vcom_ref = cvcom(q0, v0)

    # tuning variables

    # damping of the foot contact constraint
    d_bougemarte = 2e2 

    # k gain of the linear momentum error
    k_lin = 1e0
    # k gain of the qdd base error
    k_qdd_b = 1e1
    # k gain of the qdd joint error
    k_qdd_j = 1e0
    # d gain of the qdd base error
    d_qdd_b = 1e2
    # d gain of the qdd joint error
    d_qdd_j = 1e-2
    # d gain of the linear momentum eror
    d_lin = 1e-2

    # Note (SoS = sum of squares)
    # weight of each SoS body term 
    w_body = 1e0
    # weight of each SoS joint term
    w_joints = 1e0
    # weight of each SoS angular momentum term
    w_angmom = 1e0
    # weight of each SoS linear momentum term
    w_linmom = 1e0

    # weight of the body tracking objective
    beta_body = 5e2
    # weight of the joint tracking objective
    beta_joints = 1e0
    # weight of the angular momentum tracking objective
    beta_angmom = 1e0
    # weight of the linear momentum tracking objetive
    beta_linmom = 1e0
        
    i = 0
    unpause_physics_client(EmptyRequest())
    marray = Float64MultiArray()
    marray.data = np.zeros((12, 1))
    tau_publisher.publish(marray)
    while True:
        # for i in range(10):
        ms = get_model_state_client(
            GetModelStateRequest("huron_description", "ground_plane")
        )
        p = ms.pose.position
        o = ms.pose.orientation
        vl = ms.twist.linear
        al = ms.twist.angular
        q = np.vstack(
            (np.reshape(np.array([p.x, p.y, p.z, o.x, o.y, o.z, o.w]), (7, 1)), xc)
        )
        v = np.vstack(
            (np.reshape(np.array([vl.x, vl.y, vl.z, al.x, al.y, al.z]), (6, 1)), vc)
        )

        pin.forwardKinematics(model, data, q, v)
        differ = pin.difference(model, q0, q)
        qdd_b_des = k_qdd_b * np.reshape(differ[:6], (6, 1)) + d_qdd_b * np.reshape(v0[:6] - v[:6], (6,1))
        qdd_j_des = k_qdd_j * np.reshape(q0[7:] - q[7:], (12, 1)) + d_qdd_j * np.reshape(v0[6:] - v[6:], (12, 1))
        Ldot_des = k_lin * mass * (com_ref - ccom(q)) - d_lin * mass * (vcom_ref - cvcom(q, v))
        print(com_ref)
        print(ccom(q))
        Hdot_des = np.zeros((3, 1))

        opti = ca.Opti("conic")
        var_ddq = opti.variable(nv)
        var_tau = opti.variable(12)
        var_f_left = opti.variable(6)
        var_f_right = opti.variable(6)

        objective = 0
        objective = objective + (
            beta_body * ca.sumsqr(w_body * (qdd_b_des - var_ddq[:6]))
            + beta_joints * ca.sumsqr(w_joints * (qdd_j_des - var_ddq[6:]))
            + beta_angmom * ca.sumsqr(w_angmom * (Hdot_des - Hdot_fun(q, v, var_ddq)))
            + beta_linmom * ca.sumsqr(w_linmom * (Ldot_des - Ldot_fun(q, v, var_ddq)))
            + ca.sumsqr(var_tau)
            # + ca.sumsqr(var_f_left)
            # + ca.sumsqr(var_f_right)
        )
        opti.subject_to(
            acc1_constraint(q, v, var_ddq)
            == -2.0
            * np.sqrt(d_bougemarte)
            * pin.getFrameVelocity(model, data, l_ID, pin.LOCAL).vector
        )
        opti.subject_to(
            acc2_constraint(q, v, var_ddq)
            == -2.0
            * np.sqrt(d_bougemarte)
            * pin.getFrameVelocity(model, data, r_ID, pin.LOCAL).vector
        )
        opti.subject_to(
            cdynamics(
                q,
                v,
                var_ddq,
                ca.vertcat(np.zeros((6, 1)), var_tau),
                var_f_left,
                var_f_right,
            )
            == 0
        )
        opti.subject_to(var_f_left[5] >= 0)
        opti.subject_to(var_f_right[5] >= 0)
        # opti.subject_to(fric_constraint(var_f_left, var_f_right) < np.zeros((4*num_contacts,1)))
        # opti.subject_to(fric_constraint(var_f_left, var_f_right) > -np.inf*np.ones((4*num_contacts,1)))
        f_guess = ca.vertcat(0, 0, 0, 0, 0, data.mass[0] * 9.81 / 2.0)
        tau_guess = getStaticTorque(q, v, f_guess, f_guess)
        opti.set_initial(var_f_left, f_guess)
        opti.set_initial(var_f_right, f_guess)
        opti.set_initial(var_tau, tau_guess[6:])
        opti.bounded(-300, var_tau, 300)
        if i > 0:
            opti.set_initial(opti.x, prev_x)
        opti.minimize(objective)
        opti.solver("osqp", {"verbose": True})
        sol = opti.solve_limited()
        prev_x = opti.value(opti.x)
        optimized_tau = opti.value(var_tau)
        optimized_fleft = opti.value(var_f_left)
        optimized_fright = opti.value(var_f_right)
        marray = Float64MultiArray()
        marray.data = optimized_tau
        tau_publisher.publish(marray)
        print(optimized_tau)
        # print(optimized_fleft)
        # print(optimized_fright)
        # print("mass:\n",data.mass[0])
        # print("M:\n", getM(q))
        # print("C:\n", getC(q, v))
        # print("g:\n", getG(q))

    pause_physics_client(EmptyRequest())


if __name__ == "__main__":
    try:
        main()
    except rospy.ROSInterruptException:
        pass
