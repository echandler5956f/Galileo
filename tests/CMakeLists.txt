# This file is located in the tests/ directory.

# Add the test executable
add_executable(test_runner
    # List all test files here
    common/test_linalg.cpp
    common/test_polynomial.cpp
    multibody/states/test_multibody_state.cpp
    multibody/residuals/test_frame_velocity.cpp
    # Add new test files here as they are created
)

# Link the test runner against our library and Catch2
target_link_libraries(test_runner PRIVATE galileo Catch2::Catch2WithMain)

# Add the test to ctest
add_test(NAME galileo_tests COMMAND test_runner) 