cmake_minimum_required(VERSION 3.13)

project(GalileoExamples)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_custom_target(examples)

find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)
find_package(mujoco REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Threads REQUIRED)
find_package(qpOASES REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if (CMAKE_COMPILER_IS_GNUCXX)
   set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
   add_definitions ( "-DNDEBUG" )
endif (CMAKE_COMPILER_IS_GNUCXX)

include_directories(${EIGEN3_INCLUDE_DIRS} ${mujoco_INCLUDE_DIRS} ${qpOASES_INCLUDE_DIRS})

add_executable(mujoco-test mujoco-test.cpp)
target_link_libraries(mujoco-test Eigen3::Eigen ${PINOCCHIO_LIBRARIES} mujoco glfw Threads::Threads qpOASES ${CMAKE_DL_LIBS})