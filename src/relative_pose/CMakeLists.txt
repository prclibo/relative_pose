# add_library(pc_4pra pc_4pra.cpp)
# target_link_libraries(pc_4pra ${OpenCV_LIBS})
# 
# add_library(pc_5p_lih pc_5p_lih.cpp)
# target_link_libraries(pc_5p_lih lih_5pt ${OpenCV_LIBS})
# 
# add_library(pc_4pst0_nulle pc_4pst0_nulle.cpp)
# target_link_libraries(pc_4pst0_nulle ${OpenCV_LIBS})

add_library(relative_pose SHARED
    pc_4pra.cpp
    pc_5p_lih.cpp
    pc_4pst0_nulle.cpp pc_4pst0_nulle_poly.cpp
    pc_3prast0_t2d.cpp pc_3prast0_t2d_poly.cpp)
target_link_libraries(relative_pose
    lih_5pt ${OpenCV_LIBS})
