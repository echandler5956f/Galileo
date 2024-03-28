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


### Building galileo_ros for the provided examples

