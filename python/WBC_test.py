# #!/home/quant/ros_ws/src/huron_centroidal/.venv/bin/python3

# import rospy
# from std_msgs.msg import Float64MultiArray
# from sensor_msgs.msg import JointState
# from std_srvs.srv import Empty, EmptyRequest
# from gazebo_msgs.srv import GetModelState, GetModelStateRequest

# from types import SimpleNamespace

# import time
# import casadi
# import numpy as np
# import pinocchio as pin
# from pinocchio import casadi as cpin
# from pinocchio.robot_wrapper import RobotWrapper
# import numpy as np
# import casadi as ca

# import rospy
# from gazebo_msgs.srv import SpawnModel

# xc = None
# vc = None


# def joint_callback(data):
#     global xc, vc, tauc
#     xc = data.position
#     vc = data.velocity
#     # print('Got joint data')


# class SimpleInterface:
#     def __init__(self) -> None:
#         # --- Load robot model
#         builder = RobotWrapper.BuildFromURDF
#         self.robot = builder(
#             "/home/quant/ros_ws/src/HURON-Model/huron_description/urdf/huron_cheat.urdf",
#             ["/home/quant/ros_ws/src/HURON-Model/huron_description/urdf"],
#             None,
#         )
#         self.robot.q0 = np.array(
#             [
#                 0,
#                 0,
#                 1.0627,
#                 0,
#                 0,
#                 0,
#                 1,
#                 0.0000,
#                 0.0000,
#                 -0.3207,
#                 0.7572,
#                 -0.4365,
#                 0.0000,
#                 0.0000,
#                 0.0000,
#                 -0.3207,
#                 0.7572,
#                 -0.4365,
#                 0.0000,
#             ]
#         )

#         self.model = self.robot.model
#         self.data = self.model.createData()
#         pin.framesForwardKinematics(self.model, self.data, self.robot.q0)
#         pin.computeTotalMass(self.model, self.data)
#         pin.updateFramePlacements(self.model, self.data)
#         self.mass = self.data.mass[0]

#     def get_com(self, q):
#         pin.centerOfMass(self.model, self.data, q, False)
#         return self.data.com[0]

#     def get_vcom(self, q, v):
#         pin.centerOfMass(self.model, self.data, q, v, False)
#         return self.data.vcom[0]

#     def get_A(self, q):
#         pin.computeCentroidalMap(self.model, self.data, q)
#         return self.data.Ag

#     def get_Adot(self, q, v):
#         pin.computeCentroidalMapTimeVariation(self.model, self.data, q, v)
#         return self.data.dAg

#     def get_h(self, q, v):
#         return np.matmul(self.get_A(q), v)

#     def get_hdot(self, q, v, a):
#         return np.matmul(self.get_A(q), a) + np.matmul(self.get_Adot(q, v), v)

#     def get_left_foot_pos(self, q):
#         pin.forwardKinematics(self.model, self.data, q)
#         pin.updateFramePlacements(self.model, self.data)
#         return self.data.oMf[self.model.getFrameId("l_foot_v_ft_link")].translation

#     def get_right_foot_pos(self, q):
#         pin.forwardKinematics(self.model, self.data, q)
#         pin.updateFramePlacements(self.model, self.data)
#         return self.data.oMf[self.model.getFrameId("r_foot_v_ft_link")].translation

#     def compute_dynamic_components(self, q, qdot):
#         # Mass matrix
#         pin.crba(self.model, self.data, q)
#         pin.nonLinearEffects(self.model, self.data, q, qdot)
#         M_symbolic = self.data.M
#         # mass_matrix_func = ca.Function("mass_matrix_func", [q], [M_symbolic])

#         # Coriolis matrix
#         C_symbolic = ca.mtimes(self.data.C, qdot)
#         # coriolis_matrix_func = ca.Function("coriolis_matrix_func", [q, qdot], [C_symbolic])

#         # Gravity vector
#         g_sym = self.data.g
#         # gravity_vector_func = ca.Function("gravity_vector", [q], [g_sym])

#         # Torque selection matrix
#         S = np.zeros((self.model.nv, self.model.nv))
#         # 12 x 18
#         S[6:, 6:] = np.eye(self.model.nv - 6)

#         return M_symbolic, C_symbolic, g_sym, S


# def main():
#     global xc, vc, tauc

#     print("Start SimpleInterface setup")

#     q0 = np.array(
#         [
#             0,
#             0,
#             1.0627,
#             0,
#             0,
#             0,
#             1,
#             0.0000,
#             0.0000,
#             -0.3207,
#             0.7572,
#             -0.4365,
#             0.0000,
#             0.0000,
#             0.0000,
#             -0.3207,
#             0.7572,
#             -0.4365,
#             0.0000,
#         ]
#     )
#     ibrahim = SimpleInterface()

#     model = ibrahim.model
#     data = ibrahim.data
#     mass = ibrahim.mass

#     print("SimpleInterface setup done")

#     rospy.wait_for_service("/gazebo/get_model_state")
#     unpause_physics_client = rospy.ServiceProxy("/gazebo/unpause_physics", Empty)
#     pause_physics_client = rospy.ServiceProxy('/gazebo/pause_physics', Empty)

#     get_model_state_client = rospy.ServiceProxy(
#         "/gazebo/get_model_state", GetModelState
#     )
#     rospy.wait_for_service("/gazebo/spawn_urdf_model")
#     spawn_model_service = rospy.ServiceProxy("/gazebo/spawn_urdf_model", SpawnModel)

#     rospy.Subscriber("/huron/joint_states", JointState, joint_callback, queue_size=100)
#     tau_publisher = rospy.Publisher(
#         "/huron/joint_group_effort_controller/command",
#         Float64MultiArray,
#         queue_size=100,
#     )
#     rospy.init_node("ibrahim_test", anonymous=True)

#     DT = 0.08
#     contacts = [
#         SimpleNamespace(name="l_foot_v_ft_link", type=pin.ContactType.CONTACT_6D),
#         SimpleNamespace(name="r_foot_v_ft_link", type=pin.ContactType.CONTACT_6D),
#     ]

#     # initialize symbolic model and data
#     cmodel = cpin.Model(model)
#     cdata = cmodel.createData()

#     nq = model.nq
#     nv = model.nv
#     nx = nq + nv
#     ndx = 2 * nv
#     cx = casadi.SX.sym("x", nx, 1)
#     cdx = casadi.SX.sym("dx", ndx, 1)
#     cq = cx[:nq]
#     cv = cx[nq:]
#     caq = casadi.SX.sym("a", nv, 1)
#     ctauq = casadi.SX.sym("tau", nv, 1)

#     # Error in contact position
#     dpcontacts = {}
#     # Error in contact velocity
#     vcontacts = {}

#     for c in contacts:
#         c.id = model.getFrameId(c.name)
#         assert c.id < len(model.frames)
#         c.jid = model.frames[c.id].parentJoint
#         c.placement = model.frames[c.id].placement
#         c.model = pin.RigidConstraintModel(c.type, model, c.jid, c.placement)
#     contact_models = [c.model for c in contacts]

#     for c in contacts:
#         p0 = data.oMf[c.id]
#         dpcontacts[c.name] = ca.Function(
#             f"dpcontact_{c.name}",
#             [cx],
#             [np.zeros(6)],
#             ["x"],
#             ["6D_pos_error"],
#         )
#         vcontacts[c.name] = ca.Function(
#             f"vcontact_{c.namee}",
#             [cx],
#             [
#                 cpin.getFrameVelocity(
#                     cmodel, cdata, c.id, pin.LOCAL
#                 ).vector
#             ],
#             [""],
#             ["6D_vel_error"],
#         )

#     unpause_physics_client(EmptyRequest())

#     r = rospy.Rate(1 / DT)
#     i = 0
#     # for j in range(2):
#     while True:
#         seconds = rospy.get_time()
#         ms = get_model_state_client(
#             GetModelStateRequest("huron_description", "ground_plane")
#         )
#         p = ms.pose.position
#         o = ms.pose.orientation
#         vl = ms.twist.linear
#         al = ms.twist.angular
#         q = np.hstack(([p.x, p.y, p.z, o.x, o.y, o.z, o.w], xc))
#         v = np.hstack(([vl.x, vl.y, vl.z, al.x, al.y, al.z], vc))

#         # current joint angles
#         q_cur = q
#         qdot_cur = v

#         print("Start opti")

#         opti = casadi.Opti()
        
#         # simple solver options
#         opti.solver(
#             "ipopt", {"expand": False}, {"max_iter": 10}
#         )
#         if i == 0:
#             pause_physics_client(EmptyRequest())
#         # solve the problem
#         sol = opti.solve_limited()
#         # get the solution
#         # sol_us = [opti.value(var_u) for var_u in var_us]
#         # prev_lam_g = opti.value(opti.lam_g)
#         # prev_x = opti.value(opti.x)
#         # optimized_tau = opti.value(sol_us[0])

#         last_12_torques = optimized_tau

#         marray = Float64MultiArray()
#         marray.data = last_12_torques
#         tau_publisher.publish(marray)
#         print(optimized_tau)
#         i = i + 1
#         DT = rospy.get_time() - seconds


# if __name__ == "__main__":
#     try:
#         main()
#     except rospy.ROSInterruptException:
#         pass
