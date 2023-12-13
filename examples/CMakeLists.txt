cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 11)

project (examples)

add_executable(simple simple_traj_test.cpp)
add_executable(test test.cpp)
add_executable(traj_test traj_test.cpp)

target_link_libraries(simple
    PUBLIC
    galileo
)

target_link_libraries(test
    PUBLIC
    galileo
)

target_link_libraries(traj_test 
    PUBLIC
    galileo
)