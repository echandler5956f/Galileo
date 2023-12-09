# Source/Include files
include_directories(".")
include_directories("include/" ) # acrobatics includes

include_directories(${CMAKE_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
file(GLOB_RECURSE sources "src/*.cpp") # acrobatics library cpp files

# Library
add_library(galileo SHARED ${sources}) # produce a library used by interface
target_link_libraries(galileo model)
