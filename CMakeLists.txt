cmake_minimum_required(VERSION 3.13)
project(galileo)
set(PROJECT_NAME galileo)
set(PROJECT_DESCRIPTION
    "A general and efficient optimal control solver library for switched systems.")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wcast-align -Wpedantic
                    -Wdisabled-optimization -Wno-sign-conversion
                    -Wformat=2 -Winit-self -Wmissing-include-dirs 
                    -Woverloaded-virtual -Wshadow -Wsign-promo)

option(BUILD_BENCHMARK "Build the benchmarks" ON)
option(BUILD_EXAMPLES "Build the examples" ON)

option(ENABLE_VECTORIZATION_SUOPPORT 
    "Enable vectorization and further processor-related optimizations" OFF)
option(BUILD_WITH_OPENMP_SUPPORT "Build the library with OpenMP multithreading support" ON)
option(BUILD_WITH_CASADI_SUPPORT "Build the library with CasADi for AD and code generation" ON)
option(BUILD_WITH_CPPADCODEGEN_SUPPORT "Build the library with CppADCodeGen for AD and code generation" ON)
option(BUILD_WITH_COLLISION_SUPPORT 
    "Build the library with the collision support (required HPP-FCL)" ON)
option(BUILD_WITH_PINOCCHIO_SUPPORT 
    "Build the library with Pinocchio for kinematics and dynamics" ON)

if(ENABLE_VECTORIZATION_SUOPPORT)
    add_definitions(-march=native)
endif()
if(BUILD_WITH_OPENMP_SUPPORT)
  add_definitions(-D_GLIBCXX_PARALLEL)
  add_compile_options(-fopenmp -lgomp)
  find_package(OpenMP REQUIRED)
endif()
if(BUILD_WITH_CASADI_SUPPORT)
  find_package(casadi REQUIRED)
endif()
if(BUILD_WITH_CPPADCODEGEN_SUPPORT)
  find_package(cppadcg REQUIRED)
endif()
if(BUILD_WITH_COLLISION_SUPPORT)
  find_package(hpp-fcl REQUIRED)
endif()
if(BUILD_WITH_PINOCCHIO_SUPPORT)
  find_package(pinocchio REQUIRED)
endif()

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS iostreams)

add_library(galileo SHARED galileo.cc)
set_target_properties(galileo PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(third-party)
add_subdirectory(galileo)

target_link_libraries(galileo
    PUBLIC
    Eigen3::Eigen
    ${Boost_LIBRARIES})

if(BUILD_WITH_OPENMP_SUPPORT)
  target_link_libraries(galileo
    PUBLIC
    OpenMP::OpenMP_CXX)
endif()
if(BUILD_WITH_CASADI_SUPPORT)
  target_link_libraries(galileo
  PUBLIC
  casadi)
endif()
if(BUILD_WITH_CPPADCODEGEN_SUPPORT)
  target_link_libraries(galileo
  PUBLIC
  cppadcg)
endif()
if(BUILD_WITH_COLLISION_SUPPORT)
  target_link_libraries(galileo
  PUBLIC
  hpp-fcl::hpp-fcl)
endif()
if(BUILD_WITH_PINOCCHIO_SUPPORT)
  target_link_libraries(galileo
  PUBLIC
  pinocchio::pinocchio)
endif()

target_include_directories(galileo
    PUBLIC
    ${Boost_INCLUDE_DIRS})

# Install the "galileoConfig.cmake" file and export the target
install(TARGETS galileo
  EXPORT galileoTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

# Generate and install the export targets for "galileo"
install(EXPORT galileoTargets
  FILE galileoTargets.cmake
  NAMESPACE galileo::
  DESTINATION lib/cmake/galileo)

# Generate and install galileoConfig.cmake
include(CMakePackageConfigHelpers)

# Generate ${PROJECT_NAME}Config.cmake
configure_package_config_file(
  galileoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/galileoConfig.cmake
  INSTALL_DESTINATION lib/cmake/galileo)

# Install the generated "galileoConfig.cmake" file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/galileoConfig.cmake
    DESTINATION lib/cmake/galileo)
install(FILES ${CMAKE_MODULE_PATH}/FindCASADI.cmake
    DESTINATION lib/cmake/galileo/Modules)

# Uninstall target
# Source: https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

add_subdirectory(examples)
