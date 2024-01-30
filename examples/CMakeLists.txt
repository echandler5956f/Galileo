cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (examples)

add_definitions(-D_GLIBCXX_PARALLEL)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -lgomp")

add_executable(simple_test src/simple_test.cpp)
add_executable(model_building_test src/model_building_test.cpp)
add_executable(traj_test src/traj_test.cpp)

# Find the OpenMP package
find_package(OpenMP)

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

if (OpenMP_CXX_FOUND)
    # Link your target with the OpenMP library
    target_link_libraries(traj_test PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(simple_test
    PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(traj_test
    PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)