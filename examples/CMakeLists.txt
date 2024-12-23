add_executable(huron_test src/huron_test.cpp)
add_executable(go1_test src/go1_test.cpp)
add_executable(atlas_test src/atlas_test.cpp)
add_executable(simple_test src/simple_test.cpp)
add_executable(model_building_test src/model_building_test.cpp)

target_link_libraries(huron_test 
    PUBLIC
    galileo)
target_link_libraries(go1_test 
    PUBLIC
    galileo)
target_link_libraries(atlas_test 
    PUBLIC
    galileo)
target_link_libraries(simple_test
    PUBLIC
    galileo)
target_link_libraries(model_building_test
    PUBLIC
    galileo)