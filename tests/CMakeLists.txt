# Find all test source files in the tests directory recursively
file(GLOB_RECURSE TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# Add the test executable
add_executable(test_runner
    ${TEST_SOURCES}
)

# Add tests directory to include path (for test utilities)
target_include_directories(test_runner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link the test runner against our library and Catch2
target_link_libraries(test_runner PRIVATE galileo Catch2::Catch2WithMain)

# Add the test to ctest
add_test(NAME galileo_tests COMMAND test_runner)
