if(BUILD_WITH_MUJOCO_SUPPORT)
  add_executable(mujoco_test mujoco_test.cpp)
  target_link_libraries(mujoco_test PRIVATE galileo)
endif()

if(BUILD_WITH_CANDLEWICK_SUPPORT)
  find_package(example-robot-data REQUIRED)
  add_executable(candlewick_test candlewick_test.cpp)
  target_link_libraries(candlewick_test PRIVATE galileo example-robot-data::example-robot-data)
endif()

add_executable(go1_trot go1_trot.cpp)
target_link_libraries(go1_trot PRIVATE galileo)
