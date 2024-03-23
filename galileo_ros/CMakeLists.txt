cmake_minimum_required(VERSION 3.0.2)
project(galileo_ros)

add_compile_options(-std=c++17)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  sensor_msgs
  message_generation
)

## System dependencies are found with CMake's conventions
find_package(casadi REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)

## Generate messages in the 'msg' folder
add_message_files(
  FILES
  ContactPhase.msg
  ContactSequence.msg
  EnvironmentSurface.msg
  GalileoCommand.msg
  ParameterFileLocation.msg
  RobotModel.msg
)

add_service_files(
  FILES
  SolutionRequest.srv
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
 CATKIN_DEPENDS roscpp std_msgs message_runtime
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
)

## Declare a C++ library
add_library(${PROJECT_NAME}
  src/GalileoLeggedRos.cpp
)

## Declare a C++ executable
add_executable(${PROJECT_NAME}_GalileoLeggedRos 
  src/GalileoLeggedRos.cpp
)

target_link_libraries(${PROJECT_NAME}_GalileoLeggedRos  ${CASADI_LIBRARIES} Eigen3::Eigen galileo)

target_link_libraries(${PROJECT_NAME}_GalileoLeggedRos
   ${catkin_LIBRARIES}
)

add_executable(${PROJECT_NAME}_legged_test_node
  src/Nodes/galileo_legged_test_node.cpp
)

target_link_libraries(${PROJECT_NAME}_legged_test_node
   ${catkin_LIBRARIES}
)

