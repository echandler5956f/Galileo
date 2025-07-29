if(BUILD_WITH_MUJOCO_SUPPORT)
  add_executable(mujoco_test mujoco_test.cpp)
  target_link_libraries(mujoco_test PRIVATE galileo)
endif()

add_executable(go1_trot go1_trot.cpp)
target_link_libraries(go1_trot PRIVATE galileo)
