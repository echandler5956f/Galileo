cmake_minimum_required(VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (examples)

option(BUILD_SIMPLE_TESTS "Build the simple test examples" OFF)

add_definitions(-D_GLIBCXX_PARALLEL)

set(CMAKE_CXX_FLAGS "-g -O3 -fopenmp -lgomp")

# Find the OpenMP package
find_package(OpenMP)

add_executable(huron_test src/huron_test.cpp)
add_executable(go1_test src/go1_test.cpp)

target_link_libraries(huron_test 
    PUBLIC
    galileo
)
target_link_libraries(go1_test 
    PUBLIC
    galileo
)

if (OpenMP_CXX_FOUND)
    # Link your target with the OpenMP library
    target_link_libraries(huron_test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(go1_test PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(huron_test
    PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)
target_include_directories(go1_test
    PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)

if (BUILD_SIMPLE_TESTS)
    add_executable(simple_test src/simple_test.cpp)
    add_executable(model_building_test src/model_building_test.cpp)
    target_link_libraries(simple_test
        PUBLIC
        galileo
    )
    target_link_libraries(model_building_test
            PUBLIC
            galileo
        )
    target_include_directories(simple_test
        PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)
endif()