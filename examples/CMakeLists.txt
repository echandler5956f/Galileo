if(BUILD_WITH_MUJOCO_SUPPORT)
  add_executable(mujoco_test mujoco_test.cpp)
  target_link_libraries(mujoco_test PRIVATE galileo)
endif()

if(BUILD_WITH_CANDLEWICK_SUPPORT)
  find_package(example-robot-data REQUIRED)
  # add_executable(candlewick_test candlewick_test.cpp)
  # target_link_libraries(candlewick_test PRIVATE galileo example-robot-data::example-robot-data)

  # add_executable(candlewick_advanced_test candlewick_advanced_test.cpp)
  # target_link_libraries(candlewick_advanced_test PRIVATE galileo example-robot-data::example-robot-data)

  # add_executable(forward_pass_test forward_pass_test.cpp)
  # target_link_libraries(forward_pass_test PRIVATE galileo example-robot-data::example-robot-data)
endif()

# add_executable(test1 test1.cpp)
# target_link_libraries(test1 PRIVATE galileo)

# add_executable(go1_trot go1_trot.cpp)
# target_link_libraries(go1_trot PRIVATE galileo)

# add_executable(architecture_test architecture_test.cpp)
# target_link_libraries(architecture_test PRIVATE galileo)

# add_executable(pinocchio_size_check pinocchio_size_check.cpp)
# target_link_libraries(pinocchio_size_check PRIVATE galileo)

add_executable(dimension_extraction dimension_extraction.cpp)
target_link_libraries(dimension_extraction PRIVATE galileo)
