cmake_minimum_required(VERSION 3.0.2)
project(galileo_ros)

add_compile_options(-std=c++17)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  std_srvs
  message_generation
)

## System dependencies are found with CMake's conventions
find_package(casadi REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)

## Generate messages in the 'msg' folder
add_message_files(
  FILES
  EnvironmentSurface.msg
  ContactSequence.msg
  GalileoCommand.msg
  ParameterFileLocation.msg
  RobotModel.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS roscpp std_msgs std_srvs message_runtime
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
  ${GALILEO_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/GalileoLeggedRos.cpp
)

## Declare a C++ executable
add_executable(GalileoLeggedRos
  src/GalileoLeggedRos.cpp
)

target_link_libraries(GalileoLeggedRos casadi Eigen3::Eigen galileo)

target_link_libraries(GalileoLeggedRos
   ${catkin_LIBRARIES}
)

add_executable(galileo_legged_test_node
  Nodes/galileo_legged_test_node.cpp
)

target_link_libraries(galileo_legged_test_node
   ${catkin_LIBRARIES}
)