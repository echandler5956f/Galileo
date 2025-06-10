cmake_minimum_required(VERSION 3.13)

project(galileo)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(Galileo_SOURCE_DIR ${Galileo_SOURCE_DIR})
set(Galileo_BINARY_DIR ${Galileo_BINARY_DIR})

# automatically parse the version number
file(READ "${PROJECT_SOURCE_DIR}/include/galileo/macros.hpp" _galileo_version_header)
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

option(ENABLE_VECTORIZATION_SUOPPORT 
    "Enable vectorization and further processor-related optimizations" OFF)
option(BUILD_WITH_OPENMP_SUPPORT "Build the library with OpenMP multithreading support" OFF)
option(BUILD_WITH_CASADI_SUPPORT "Build the library with CasADi for AD and code generation" OFF)
option(BUILD_WITH_CPPADCODEGEN_SUPPORT "Build the library with CppADCodeGen for AD and code generation" OFF)

add_compile_options(-Wall -Wextra -Wcast-align -Wpedantic
                    -Wdisabled-optimization -Wno-sign-conversion
                    -Wformat=2 -Winit-self -Wmissing-include-dirs 
                    -Woverloaded-virtual -Wshadow -Wsign-promo)

if(ENABLE_VECTORIZATION_SUOPPORT)
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
find_package(glfw3 REQUIRED)
find_package(Threads REQUIRED)
find_package(qpOASES REQUIRED)

include_directories(${GALILEO3_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${mujoco_INCLUDE_DIRS} ${qpOASES_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

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

add_subdirectory(examples EXCLUDE_FROM_ALL)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
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