cmake_minimum_required(VERSION 3.13)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_COLOR_MAKEFILE ON)
project(Galileo)

find_package(casadi)
find_package(pinocchio)
find_package (Eigen3 REQUIRED NO_MODULE)

# libgalileo
add_library(galileo SHARED galileo.cc)
set_target_properties(galileo PROPERTIES LINKER_LANGUAGE CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb -Wall \
    -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2 \
    -Winit-self -Wmissing-include-dirs -Woverloaded-virtual \
    -Wshadow -Wsign-promo")

set(CMAKE_C_FLAGS "-O3 -ggdb -std=gnu99 -I.")

add_subdirectory(galileo)

target_link_libraries(galileo
    PUBLIC
    casadi
    pinocchio
    Eigen3::Eigen
)

# Install the "galileoConfig.cmake" file and export the target
install(TARGETS galileo
  EXPORT galileoTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  # FILE_SET HEADERS DESTINATION include
)

# Generate and install the export targets for "galileo"
install(EXPORT galileoTargets
  FILE galileoTargets.cmake
  NAMESPACE galileo::
  DESTINATION lib/cmake/galileo
)

# Generate and install galileoConfig.cmake
include(CMakePackageConfigHelpers)

# Generate ${PROJECT_NAME}Config.cmake
configure_package_config_file(
    galileoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/galileoConfig.cmake
  INSTALL_DESTINATION lib/cmake/galileo
)
# Install the generated "galileoConfig.cmake" file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/galileoConfig.cmake
    DESTINATION lib/cmake/galileo
)
install(
    FILES ${CMAKE_MODULE_PATH}/FindCASADI.cmake
    DESTINATION lib/cmake/galileo/Modules
)

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