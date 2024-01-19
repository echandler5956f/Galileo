cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (examples)

set(CMAKE_CXX_FLAGS "-O3")

add_executable(rho_beta_epsilon src/rho_beta_epsilon.cpp)
add_executable(simple_test src/simple_test.cpp)
add_executable(model_building_test src/model_building_test.cpp)
add_executable(traj_test src/traj_test.cpp)

target_link_libraries(rho_beta_epsilon
    PUBLIC
    galileo
)

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

target_include_directories(traj_test
    PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)