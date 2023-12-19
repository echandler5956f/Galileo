cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (examples)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb -Wall \
    -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2 \
    -Winit-self -Wmissing-include-dirs -Woverloaded-virtual \
    -Wshadow -Wsign-promo")

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