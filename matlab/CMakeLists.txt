set(MEX_INCLUDES
    "-I$<JOIN:$<TARGET_PROPERTY:relative_pose,INCLUDE_DIRECTORIES>, -I>")
set(MEX_LINKS
    "-l$<JOIN:$<TARGET_PROPERTY:relative_pose,LINK_LIBRARIES>, -l>")
set(MEX_LD
    "-L$<TARGET_FILE_DIR:lih_5pt> -L$<TARGET_FILE_DIR:relative_pose> -L$<JOIN:${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}, -L>")
set(SOURCE_FILE
    ${CMAKE_CURRENT_SOURCE_DIR}/wrappers.cpp)
set(MATLAB_BINARY
    ${Matlab_MAIN_PROGRAM} -nodisplay -nosplash -nodesktop)

file(GLOB WRAPPER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(WRAPPER_FILE ${WRAPPER_FILES})
    get_filename_component(WRAPPER_NAME ${WRAPPER_FILE} NAME_WE)
    set(COMPILE_SCRIPT "${COMPILE_SCRIPT}mexOpenCV ${WRAPPER_FILE} ${MEX_INCLUDES} ${MEX_LD} ${MEX_LINKS} -lrelative_pose,\n")
endforeach(WRAPPER_FILE)

file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mexOpenCV_compile.m
    CONTENT "${COMPILE_SCRIPT}")

add_custom_target(matlab
    COMMAND ${MATLAB_BINARY} -r
    "try, run('${CMAKE_CURRENT_BINARY_DIR}/mexOpenCV_compile.m'), catch ME, disp(ME.message), exit(1), end, exit(0);"
    VERBATIM
    DEPENDS relative_pose)

# https://cmake.org/cmake/help/v3.3/module/FindMatlab.html#command:matlab_add_mex
# https://cmake.org/cmake/help/v3.3/manual/cmake-generator-expressions.7.html
# https://cmake.org/cmake/help/v3.0/variable/CMAKE_LANG_IMPLICIT_LINK_DIRECTORIES.html

# set(mexOpenCV_COMMAND
#     "${Matlab_MAIN_PROGRAM} -nodisplay -nosplash -nodesktop -r 'mexOpenCV ${CMAKE_CURRENT_SOURCE_DIR}/wrappers.cpp,exit'")
# 
# https://stackoverflow.com/questions/32183975/how-to-print-all-the-properties-of-a-target-in-cmake
# get_property(INCLUDE_DIRS TARGET relative_pose PROPERTY INCLUDE_DIRECTORIES)
# get_property(LINK_LIBS TARGET relative_pose PROPERTY LINK_LIBRARIES)
# 
# message(STATUS "===========${mexOpenCV_COMMAND} ${LINK_LIBS} ")
# set(mexOpenCV_COMMAND
#     "${MATLAB_BINARY} -r 'mexOpenCV -v ${SOURCE_FILE} ${MEX_INCLUDES} ${MEX_LD} ${MEX_LINKS} -lrelative_pose,exit'")


