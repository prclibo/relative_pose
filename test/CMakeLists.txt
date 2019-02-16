include_directories(../third_party/googletest-release-1.8.1/googletest/include)

add_executable(test_all test_all.cpp)
target_link_libraries(test_all relative_pose ${OpenCV_LIBS})

add_executable(test_pc_4pra test_pc_4pra.cpp)
target_link_libraries(test_pc_4pra
    relative_pose gtest_all gtest_main ${OpenCV_LIBS})
