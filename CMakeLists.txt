cmake_minimum_required (VERSION 2.8.11)
project(relative_pose)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-O2)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include)
add_subdirectory(third_party)

add_subdirectory(src/relative_pose)
add_subdirectory(test)

find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MX_LIBRARY)
include_directories(${Matlab_INCLUDE_DIRS})
add_subdirectory(matlab)
message(STATUS "${Matlab_LIBRARIES}")
message(STATUS "${Matlab_MX_LIBRARY}")
