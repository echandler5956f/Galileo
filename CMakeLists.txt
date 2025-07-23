cmake_minimum_required(VERSION 3.13)

project(galileo)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(Galileo_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
set(Galileo_BINARY_DIR "${CMAKE_BINARY_DIR}")

# automatically parse the version number
file(READ "${PROJECT_SOURCE_DIR}/include/galileo/common/meta/macros.hpp" _galileo_version_header)
string(REGEX MATCH "define[ \t]+GALILEO_WORLD_VERSION[ \t]+([0-9]+)" _galileo_world_version_match "${_galileo_version_header}")
set(GALILEO_WORLD_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GALILEO_MAJOR_VERSION[ \t]+([0-9]+)" _galileo_major_version_match "${_galileo_version_header}")
set(GALILEO_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GALILEO_MINOR_VERSION[ \t]+([0-9]+)" _galileo_minor_version_match "${_galileo_version_header}")
set(GALILEO_MINOR_VERSION "${CMAKE_MATCH_1}")
set(GALILEO_VERSION_NUMBER ${GALILEO_WORLD_VERSION}.${GALILEO_MAJOR_VERSION}.${GALILEO_MINOR_VERSION})

include(CheckCXXCompilerFlag)
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

option(ENABLE_VECTORIZATION_SUPPORT
    "Enable vectorization and further processor-related optimizations" OFF)
option(BUILD_WITH_OPENMP_SUPPORT "Build the library with OpenMP multithreading support" OFF)
option(GALILEO_BUILD_TESTS "Build the tests for galileo" ON)
option(GALILEO_BUILD_EXAMPLES "Build the examples for galileo" ON)

add_compile_options(-Wall -Wextra -Wcast-align -Wdisabled-optimization
                    -Wno-sign-conversion -Wformat=2 -Winit-self
                    -Wmissing-include-dirs -Woverloaded-virtual -Wshadow
                    -Wsign-promo -Wno-unused-parameter)

if(ENABLE_VECTORIZATION_SUPPORT)
    add_definitions(-march=native)
endif()
if(BUILD_WITH_OPENMP_SUPPORT)
  add_definitions(-D_GLIBCXX_PARALLEL)
  add_compile_options(-fopenmp -lgomp)
  find_package(OpenMP REQUIRED)
endif()

find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)

find_package(mujoco REQUIRED)
find_package(GLFW3 REQUIRED)
find_package(Threads REQUIRED)

include_directories(${GALILEO_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${mujoco_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

set(INCLUDE_INSTALL_DIR
    "${CMAKE_INSTALL_INCLUDEDIR}/galileo"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where Galileo header files are installed"
    )
set(CMAKEPACKAGE_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/galileo/cmake"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where GalileoConfig.cmake is installed"
    )
set(PKGCONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/pkgconfig"
    CACHE PATH "The directory relative to CMAKE_PREFIX_PATH where galileo.pc is installed"
    )

include(RegexUtils)
test_escape_string_as_regex()

install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/galileo/"
    DESTINATION "${INCLUDE_INSTALL_DIR}"
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx" PATTERN ".*" EXCLUDE
)

# Install test resources if they exist
if(GALILEO_BUILD_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/resources")
    install(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/"
        DESTINATION "${CMAKE_INSTALL_DATADIR}/galileo/resources"
        FILES_MATCHING
            PATTERN "*.urdf"
            PATTERN "*.xacro"
            PATTERN "*.dae"
            PATTERN "*.stl"
            PATTERN "*.obj"
            PATTERN "*.mesh"
    )
endif()

# Package configuration variables
set (GALILEO_VERSION_STRING ${GALILEO_VERSION_NUMBER} )
set (GALILEO_VERSION_MAJOR ${GALILEO_WORLD_VERSION} )
set (GALILEO_VERSION_MINOR ${GALILEO_MAJOR_VERSION} )
set (GALILEO_VERSION_PATCH ${GALILEO_MINOR_VERSION} )
set (GALILEO_DEFINITIONS "")
set(GALILEO_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
set(GALILEO_ROOT_DIR ${CMAKE_INSTALL_PREFIX})

include(CMakePackageConfigHelpers)

add_library(galileo INTERFACE)

target_compile_definitions(galileo INTERFACE ${GALILEO_DEFINITIONS})
target_include_directories(galileo INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)

target_compile_features(galileo INTERFACE cxx_std_20)

target_link_libraries(galileo INTERFACE
  pinocchio_default
  pinocchio_parsers
  pinocchio_cppad
  pinocchio_cppadcg
  Eigen3::Eigen
  Threads::Threads
  mujoco
  glfw
)

set_target_properties(galileo PROPERTIES EXPORT_NAME Galileo)

install(TARGETS galileo EXPORT GalileoTargets)

configure_package_config_file (
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GalileoConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GalileoConfig.cmake
    PATH_VARS GALILEO_INCLUDE_DIR GALILEO_ROOT_DIR
    INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

set (_GALILEO_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})
unset (CMAKE_SIZEOF_VOID_P)
write_basic_package_version_file (GalileoConfigVersion.cmake
                                VERSION ${GALILEO_VERSION_NUMBER}
                                COMPATIBILITY SameMajorVersion)
set (CMAKE_SIZEOF_VOID_P ${_Galileo_CMAKE_SIZEOF_VOID_P})

export(TARGETS galileo NAMESPACE Galileo:: FILE GalileoTargets.cmake)
export(PACKAGE Galileo)
install(EXPORT GalileoTargets NAMESPACE Galileo:: DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UseGalileo.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/GalileoConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/GalileoConfigVersion.cmake
          DESTINATION ${CMAKEPACKAGE_INSTALL_DIR} )

# Add uninstall target
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GalileoUninstall.cmake)

if(GALILEO_BUILD_TESTS)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.2
  )
  FetchContent_MakeAvailable(Catch2)

    # Add the GALILEO_TESTING definition for test targets
  add_compile_definitions(GALILEO_TESTING)

  add_subdirectory(tests)
endif()

if(GALILEO_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
