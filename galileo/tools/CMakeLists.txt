FILE(GLOB_RECURSE SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

FILE(GLOB_RECURSE HEADER_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
)

FILE(GLOB_RECURSE TEST_SOURCE_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/test/*.cc
)

set(INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include/")

target_sources(galileo
        PRIVATE
        ${SOURCE_FILES}
)
target_include_directories(galileo
        PUBLIC
        $<BUILD_INTERFACE:${INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>
)

target_sources(galileo_tests
  PRIVATE
  ${TEST_SOURCE_FILES}
)
target_include_directories(galileo_tests
        PUBLIC
        $<BUILD_INTERFACE:${INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>
)

install(DIRECTORY ${INCLUDE_DIRS} DESTINATION include)