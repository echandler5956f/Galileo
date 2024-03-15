#!/home/quant/ros_ws/src/Galileo/.venv/bin/python3

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
    model_name_ = "go1"
    joint_names_ = ["LF_HAA", "LF_HFE", "LF_KFE", "RF_HAA", "RF_HFE", "RF_KFE", "LH_HAA", "LH_HFE", "LH_KFE", "RH_HAA", "RH_HFE", "RH_KFE"]
    pose_ = Pose(position = Point(x = 0, y = 0, z = 0.339), orientation = Quaternion(x = 0, y = 0, z = 0, w = 1))
    twist_ = Twist(linear = Vector3(x = 0, y = 0, z = 0), angular = Vector3(x = 0, y = 0, z = 0))
    joint_positions_ = [0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3]
    response2 = sms(SetModelStateRequest(ModelState(model_name = model_name_, pose = pose_, twist = twist_)))
    pause_physics_client(EmptyRequest())
    rospy.sleep(2)
    response1 = smc(SetModelConfigurationRequest(model_name = model_name_, joint_names = joint_names_, joint_positions = joint_positions_))
    print(response1.success)
    print(response1.status_message)
    print(response2.success)
    print(response2.status_message)

if __name__ == '__main__':
    try:
        main()
    except rospy.ROSInterruptException:
        pass