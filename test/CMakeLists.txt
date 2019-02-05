add_executable(test_all test_all.cpp)
target_link_libraries(test_all relative_pose ${OpenCV_LIBS})
