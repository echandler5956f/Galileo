cmake_minimum_required(VERSION 3.10.2)
project(test_env)

find_package(casadi)
find_package(pinocchio)
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS iostreams unit_test_framework)
find_package(OpenMP)

add_executable(casadi_demo casadi_demo.cpp)
target_link_libraries(casadi_demo casadi OpenMP::OpenMP_CXX pinocchio Eigen3::Eigen ${Boost_LIBRARIES})

include_directories(${Boost_INCLUDE_DIRS})
add_executable(casadi_demo2 casadi_demo2.cpp)
target_link_libraries(casadi_demo2 casadi OpenMP::OpenMP_CXX pinocchio Eigen3::Eigen ${Boost_LIBRARIES})

# add_executable(casadi_demo3 casadi_demo3.cpp)
# target_link_libraries(casadi_demo3 casadi pinocchio Eigen3::Eigen)

# add_executable(full_demo full_demo.cpp)
# target_link_libraries(full_demo casadi pinocchio Eigen3::Eigen Boost::unit_test_framework)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -O3")

install(TARGETS casadi_demo RUNTIME DESTINATION . LIBRARY DESTINATION .)

install(TARGETS casadi_demo2 RUNTIME DESTINATION . LIBRARY DESTINATION .)

# install(TARGETS casadi_demo3 RUNTIME DESTINATION . LIBRARY DESTINATION .)