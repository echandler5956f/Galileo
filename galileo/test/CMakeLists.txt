add_executable(traj_test traj_test.cpp)
add_executable(test test.cpp)
target_link_libraries(traj_test galileo)
target_link_libraries(test galileo)