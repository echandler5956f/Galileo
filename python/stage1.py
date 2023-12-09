#!/home/quant/ros_ws/src/huron_centroidal/.venv/bin/python3

import rospy
from std_msgs.msg import Float64MultiArray
from sensor_msgs.msg import JointState
from std_srvs.srv import Empty, EmptyRequest
from gazebo_msgs.srv import GetModelState, GetModelStateRequest

import time
import casadi
import numpy as np
import pinocchio as pin
from pinocchio import casadi as cpin
from pinocchio.robot_wrapper import RobotWrapper

xc = None
vc = None

def joint_callback(data):
    global xc, vc, tauc
    xc = data.position
    vc = data.velocity

class SimpleInterface:
    def __init__(self) -> None:
        # --- Load robot model
        builder = RobotWrapper.BuildFromURDF
        self.robot = builder(
            "resources/urdf/huron_cheat.urdf",
            ["resources"],
            None,
        )
        self.robot.q0 = np.array([0, 0, 1.0627, 0, 0, 0, 1,
                            0.0000,  0.0000, -0.3207, 0.7572, -0.4365,  0.0000,
                            0.0000,  0.0000, -0.3207, 0.7572, -0.4365,  0.0000])
        
        self.model = self.robot.model
        self.data = self.model.createData()
        pin.framesForwardKinematics(self.model, self.data, self.robot.q0)
        pin.computeTotalMass(self.model, self.data)
        pin.updateFramePlacements(self.model, self.data)

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
        return np.matmul(self.get_A(q), a) +  np.matmul(self.get_Adot(q, v), v)
    
    def get_left_foot_pos(self, q):
        pin.forwardKinematics(self.model, self.data, q)
        pin.updateFramePlacements(self.model, self.data)
        return self.data.oMf[self.model.getFrameId("l_foot_v_ft_link")].translation

    def get_right_foot_pos(self, q):
        pin.forwardKinematics(self.model, self.data, q)
        pin.updateFramePlacements(self.model, self.data)
        return self.data.oMf[self.model.getFrameId("r_foot_v_ft_link")].translation


def main():
    global xc, vc, tauc

    q0 = np.array([0, 0, 1.0627, 0, 0, 0, 1,
                    0.0000,  0.0000, -0.3207, 0.7572, -0.4365,  0.0000,
                    0.0000,  0.0000, -0.3207, 0.7572, -0.4365,  0.0000])
    ibrahim = SimpleInterface()

    rospy.wait_for_service('/gazebo/unpause_physics')
    rospy.wait_for_service('/gazebo/pause_physics')
    rospy.wait_for_service('gazebo/get_model_state')
    unpause_physics_client = rospy.ServiceProxy('/gazebo/unpause_physics', Empty)
    # pause_physics_client = rospy.ServiceProxy('/gazebo/pause_physics', Empty)
    get_model_state_client = rospy.ServiceProxy('/gazebo/get_model_state', GetModelState)

    rospy.Subscriber("/huron/joint_states", JointState, joint_callback, queue_size=100)

    rospy.init_node('ibrahim_test', anonymous=True)
    unpause_physics_client(EmptyRequest())
    time.sleep(2)
    while True:
        ms = get_model_state_client(GetModelStateRequest("huron", "ground_plane"))
        p = ms.pose.position
        o = ms.pose.orientation
        vl = ms.twist.linear
        al = ms.twist.angular
        q = np.hstack(([p.x, p.y, p.z, o.x, o.y, o.z, o.w], xc))
        v = np.hstack(([vl.x, vl.y, vl.z, al.x, al.y, al.z], vc))
        print(ibrahim.get_h(q.flatten(), v.flatten()))
        # print(ibrahim.get_hdot(q.flatten(), v.flatten()))

if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass