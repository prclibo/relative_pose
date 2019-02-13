add_library(lih_5pt
    5pt-6pt-Li-Hartley/Ematrix_5pt.cc
    5pt-6pt-Li-Hartley/Ematrix_6pt.cc
    5pt-6pt-Li-Hartley/polydet.cc
    5pt-6pt-Li-Hartley/polyquotient.cc
    5pt-6pt-Li-Hartley/sturm.cc)

include_directories(googletest-release-1.8.1/googletest/include)
include_directories(googletest-release-1.8.1/googletest/)

add_library(gtest_main googletest-release-1.8.1/googletest/src/gtest_main.cc)
add_library(gtest_all googletest-release-1.8.1/googletest/src/gtest-all.cc)

target_link_libraries(gtest_all pthread m)
