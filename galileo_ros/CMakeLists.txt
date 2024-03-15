cmake_minimum_required(VERSION 3.0.2)
project(galileo_ros)

add_compile_options(-std=c++17)

### UNITREE SDK ###

message("-- CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64.*")
  set(ARCH amd64)
endif()
if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "aarch64.*")
  set(ARCH arm64)
endif()


link_directories(lib/cpp/${ARCH})


## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  sensor_msgs
  message_generation
)

set(UNITREE_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/../unitree_legged_sdk/include)
set(UNITREE_DIR /usr/local/lib/unitree)

## System dependencies are found with CMake's conventions
find_package(casadi REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)
find_package(unitree_legged_sdk REQUIRED)

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
  ${UNITREE_INCLUDE_DIRS}
)



## Declare a C++ library
add_library(${PROJECT_NAME}
  src/galileo_legged_ros_implementation.cpp
)

## Declare a C++ executable

# legged node
add_executable(${PROJECT_NAME}_legged_node 
  src/galileo_legged_ros_implementation.cpp
)

target_link_libraries(${PROJECT_NAME}_legged_node casadi Eigen3::Eigen galileo)

target_link_libraries(${PROJECT_NAME}_legged_node
   ${catkin_LIBRARIES}
)

# publish go1 information
add_executable(${PROJECT_NAME}_go1_publisher 
  src/helper_nodes/galileo_legged_example.cpp
)

target_link_libraries(${PROJECT_NAME}_go1_publisher
   ${catkin_LIBRARIES}
)

# poll solution
add_executable(${PROJECT_NAME}_go1_poll_solution 
  src/helper_nodes/galileo_legged_poll_solution.cpp
)

target_link_libraries(${PROJECT_NAME}_go1_poll_solution 
   ${catkin_LIBRARIES}
)

# commanding controller
add_executable(${PROJECT_NAME}_commanding_controller 
  src/commanding_controller.cpp
)



## UNITREE SDK
find_package(catkin QUIET)
if(${catkin_FOUND})
  catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_SOURCE_DIR}/lib/cpp/${ARCH}/libunitree_legged_sdk.a
  )
endif()

set(EXTRA_LIBS -pthread libunitree_legged_sdk.a)



# Add  static library
add_library(unitree_legged_sdk STATIC IMPORTED)

set_property(TARGET unitree_legged_sdk PROPERTY IMPORTED_LOCATION ${UNITREE_DIR}/libunitree_legged_sdk.a)

target_link_libraries(${PROJECT_NAME}_commanding_controller casadi Eigen3::Eigen galileo unitree_legged_sdk)

target_link_libraries(${PROJECT_NAME}_commanding_controller 
   ${catkin_LIBRARIES}
   ${EXTRA_LIBS} 
)

target_include_directories(${PROJECT_NAME}_commanding_controller PRIVATE ${UNITREE_INCLUDE_DIR})

install(TARGETS ${PROJECT_NAME}_commanding_controller DESTINATION ${CMAKE_INSTALL_PREFIX})


install(DIRECTORY lib/cpp/${ARCH}/
  DESTINATION lib/unitree
  USE_SOURCE_PERMISSIONS)