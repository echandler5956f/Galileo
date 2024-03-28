# Galileo ROS
*A ROS package for using Galileo for Legged Robots*


## Running the Legged Examples

### Building galileo_ros for the provided examples
Currently, only legged robots are supported with galileo_ros, through "GalileoLeggedRos". 

Examples exist for the Unitree Go1, and WPI's HURON robot. Mesh files and solver specific params can be found in "galileo/resources/". Before running either of these examples make sure to build the corresponding resource file with catkin build like so from inside a workspace -- 


```bash
 catkin build galileo_go1_description 
 ```

or 

```bash
 catkin build galileo_huron_description 
 ```

The ROS package can be built similarly, with 


```bash
 catkin build galileo_ros
 ```

upon which you can run the launch files. 


### galileo_ros for the provided examples

```bash
 roslaunch galileo_ros <robot_name>_galileo_ros.launch
 ```
runs the solver, and  

```bash
 roslaunch galileo_ros <robot_name>_galileo_ros_rviz.launch
 ```

the visualizer. 


In the case that you want to change the problem parameters, such as target location, configuration, or the phase sequence, make sure to look at the "problem_parameters" file under resources.

<p align="center">
  <a> <b>end_effector_names:</b> The names of the end effector frames in the URDF.</a> •
  <a> <b>knot_num:</b> the "fidelity", or number of knot points, per phase. The size of "knot_num" must be the number of phases </a> •
  <a><b>knot_time:</b> the duration of the corresponding phase. The size of "knot_time" must be the number of phases </a> •
  <a><b>contact_combinations</b> The contact surface each end effector is in contact with during some phase. "-1" represents an end-effector not in contact. This is of size number_of_phases x number_of_end_effectors.  </a> •
  <a><b>q0</b> The initial configuration of the robot. q0[0:2] is the initial global position, q0[3:6] is the quaternion representing rotation, and q0[7:end] is the configuration for each joint. </a> •
  <a><b>qf</b> The goal configuration of the robot. qf[0:2] is the final global position, qf[3:6] is the quaternion representing rotation, and qf[7:end] is the configuration for each joint. </a> •
</p>
<br/>


