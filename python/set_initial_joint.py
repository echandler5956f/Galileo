#!/home/quant/ros_ws/src/huron_centroidal/.venv/bin/python3

from gazebo_msgs.srv import SetModelConfiguration, SetModelConfigurationRequest, SetModelState, SetModelStateRequest
from gazebo_msgs.msg import ModelState
from geometry_msgs.msg import Point, Quaternion, Pose, Twist, Vector3
from std_srvs.srv import Empty, EmptyRequest
import rospy

def main():
    rospy.wait_for_service('/gazebo/set_model_configuration')
    rospy.wait_for_service('/gazebo/set_model_state')
    rospy.wait_for_service('/gazebo/unpause_physics')
    rospy.wait_for_service('/gazebo/pause_physics')

    smc = rospy.ServiceProxy('/gazebo/set_model_configuration', SetModelConfiguration)
    sms = rospy.ServiceProxy('/gazebo/set_model_state', SetModelState)
    unpause_physics_client = rospy.ServiceProxy('/gazebo/unpause_physics', Empty)
    pause_physics_client = rospy.ServiceProxy('/gazebo/pause_physics', Empty)
    
    model_name_ = "huron_description"
    urdf_param_name_ = "$(find xacro)/xacro --inorder '$(find huron_description)/urdf/huron.xacro'"
    joint_names_ = ["l_hip_yaw_joint", "l_hip_roll_joint", "l_hip_pitch_joint", "l_knee_pitch_joint", "l_ankle_pitch_joint", "l_ankle_roll_joint",
                    "r_hip_yaw_joint", "r_hip_roll_joint", "r_hip_pitch_joint", "r_knee_pitch_joint", "r_ankle_pitch_joint", "r_ankle_roll_joint"]
    pose_ = Pose(position = Point(x = 0, y = 0, z = 1.1064), orientation = Quaternion(x = 0, y = 0, z = 0, w = 1))
    twist_ = Twist(linear = Vector3(x = 0, y = 0, z = 0), angular = Vector3(x = 0, y = 0, z = 0))
    joint_positions_ = [0.0, 0.0, -0.1672, 0.3922, -0.2251, 0.0,
                        0.0, 0.0, -0.1672, 0.3922, -0.2251, 0.0]
    pause_physics_client(EmptyRequest())
    response2 = sms(SetModelStateRequest(ModelState(model_name = model_name_, pose = pose_, twist = twist_)))
    rospy.sleep(2)
    response1 = smc(SetModelConfigurationRequest(model_name = model_name_, urdf_param_name = urdf_param_name_, joint_names = joint_names_, joint_positions = joint_positions_))
    # print(response.success)
    # print(response.status_message)

if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass