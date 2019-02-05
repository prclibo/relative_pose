cmake_minimum_required (VERSION 2.8.11)
project(relative_pose)

set (CMAKE_CXX_STANDARD 11)

find_package(OpenCV REQUIRED)

include_directories(include)
add_subdirectory(src/relative_pose)
add_subdirectory(test)
