find_package(casadi)
find_package(pinocchio)
find_package (Eigen3 REQUIRED NO_MODULE)

add_subdirectory(model)
add_subdirectory(variables)
add_subdirectory(test)