cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 11)

project (examples)

add_executable(simple_test simple_test.cpp)
add_executable(model_building_test model_building_test.cpp)
add_executable(traj_test traj_test.cpp)

target_link_libraries(simple_test
    PUBLIC
    galileo
)

target_link_libraries(model_building_test
    PUBLIC
    galileo
)

target_link_libraries(traj_test 
    PUBLIC
    galileo
)