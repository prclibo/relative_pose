add_library(pc_4pra pc_4pra.cpp)
target_link_libraries(pc_4pra ${OpenCV_LIBS})

add_library(pc_5p_lih pc_5p_lih.cpp)
target_link_libraries(pc_5p_lih lih_5pt ${OpenCV_LIBS})

add_library(pc_4pst0_nulle pc_4pst0_nulle.cpp)
target_link_libraries(pc_4pst0_nulle ${OpenCV_LIBS})

