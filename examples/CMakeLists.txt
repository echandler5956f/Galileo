cmake_minimum_required(VERSION 3.13)

project(GalileoExamples)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_custom_target(examples)

find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)
find_package(mujoco REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Threads REQUIRED)
find_package(qpOASES REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wcast-align -Wpedantic
                    -Wdisabled-optimization -Wno-sign-conversion
                    -Wformat=2 -Winit-self -Wmissing-include-dirs 
                    -Woverloaded-virtual -Wshadow -Wsign-promo)

include_directories(${EIGEN3_INCLUDE_DIRS} ${mujoco_INCLUDE_DIRS} ${qpOASES_INCLUDE_DIRS})

add_executable(mujoco-test mujoco-test.cpp)
target_link_libraries(mujoco-test Eigen3::Eigen ${PINOCCHIO_LIBRARIES} mujoco glfw Threads::Threads qpOASES ${CMAKE_DL_LIBS})

add_executable(main main.cpp)
target_link_libraries(main Eigen3::Eigen ${PINOCCHIO_LIBRARIES} ${CMAKE_DL_LIBS})