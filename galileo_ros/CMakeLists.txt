cmake_minimum_required(VERSION 3.0.2)
project(galileo_ros)

add_compile_options(-std=c++17)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  rospy
  roscpp
  std_msgs
  sensor_msgs
  message_generation
  controller_manager
  gazebo_ros
  urdf
)

## System dependencies are found with CMake's conventions
find_package(casadi REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)
find_package(gazebo REQUIRED)

## Generate messages in the 'msg' folder
add_message_files(
  FILES
  ModelLocation.msg
  RobotSolution.msg
  RobotCommand.msg
)

add_service_files(
  FILES
  DesiredStateInputCmd.srv
)

## Generate added messages and services with any dependencies listed here
generate_messages(
  DEPENDENCIES
  std_msgs
  sensor_msgs
)

###################################
## catkin specific configuration ##
###################################
catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS rospy roscpp std_msgs message_runtime controller_manager gazebo_ros urdf
)

###########
## Build ##
###########

## Specify additional locations of header files
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${GALILEO_INCLUDE_DIRS}
  ${GAZEBO_INCLUDE_DIRS}
)

## Declare a C++ library
add_library(${PROJECT_NAME}
  src/galileo_legged_ros_implementation.cpp
)

## Declare a C++ executable
add_executable(${PROJECT_NAME}_legged_node 
  src/galileo_legged_ros_implementation.cpp
)

target_link_libraries(${PROJECT_NAME}_legged_node casadi Eigen3::Eigen galileo)

target_link_libraries(${PROJECT_NAME}_legged_node
   ${catkin_LIBRARIES}
)

add_executable(${PROJECT_NAME}_go1_publisher 
  src/helper_nodes/galileo_legged_example.cpp
)

target_link_libraries(${PROJECT_NAME}_go1_publisher
   ${catkin_LIBRARIES}
)


add_executable(switch_controller_client
  src/helper_nodes/switch_controller_client.cpp
)

target_link_libraries(switch_controller_client
   ${catkin_LIBRARIES}
)