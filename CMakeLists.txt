cmake_minimum_required(VERSION 3.10.2)
project(Galileo)

set(CMAKE_COLOR_MAKEFILE ON)
# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# set(CMAKE_VERBOSE_MAKEFILE ON)

# CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
#     ${CMAKE_BINARY_DIR}/Configuration.h)

set(CMAKE_CXX_FLAGS "-O3 -no-pie -ggdb -Wall \
    -Wextra -Wcast-align -Wdisabled-optimization -Wformat=2 \
    -Winit-self -Wmissing-include-dirs -Woverloaded-virtual \
    -Wshadow -Wsign-promo")

set(CMAKE_C_FLAGS "-O3  -ggdb -std=gnu99 -I.")
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(galileo)