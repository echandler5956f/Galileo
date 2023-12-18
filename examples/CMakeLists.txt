cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (examples)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "-O3 -ggdb -std=gnu99 -I.")

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