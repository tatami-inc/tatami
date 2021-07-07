#########################################
# Building the char2double demonstration.

add_executable(
    char2double
    src/char2double.cpp
)

target_link_libraries(
    char2double 
    tatami
)

#########################################
# Building the colsums demonstration.

add_executable(
    colsums 
    src/colsums.cpp
)

target_link_libraries(
    colsums 
    tatami
)

#########################################
# Building the OpenMP demonstration.

add_executable(
    parallel 
    src/parallel.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(parallel PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
    parallel 
    PUBLIC
    tatami
)
