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

#########################################
# Checking sparse workspace performance.

add_executable(
    sparse_extractor 
    src/sparse_extractor.cpp
)

target_link_libraries(
    sparse_extractor
    tatami
)

#########################################
# Checking HDF5 reading performance.

set(TEST_HDF5 OFF CACHE BOOL "Test HDF5 support")
if (TEST_HDF5)
    add_executable(
        read_h5
        src/read_h5.cpp
    )

    target_link_libraries(
        read_h5
        tatami
    )

    enable_language(C)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    target_link_libraries(read_h5 hdf5::hdf5 hdf5::hdf5_cpp)

    target_compile_options(
        read_h5 PUBLIC -O2
    )
endif()
