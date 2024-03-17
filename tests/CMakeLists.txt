include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/afd902e992b720d1b3e106bc5e425a5768872265.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Avoid installing GoogleTest when installing this project.
option(INSTALL_GTEST "Enable installation of googletest." OFF)

FetchContent_MakeAvailable(googletest)

enable_testing()

# Main test executable.
add_executable(
    libtest 

    src/dense/DenseMatrix.cpp
    src/dense/convert_to_dense.cpp

    src/other/DelayedBind.cpp
    src/other/DelayedTranspose.cpp
    src/other/DelayedCast.cpp

    src/subset/DelayedSubset.cpp
    src/subset/DelayedSubsetBlock.cpp

##    src/isometric/unary/arith_vector_helpers.cpp
##    src/isometric/unary/arith_scalar_helpers.cpp
##    src/isometric/unary/math_helpers.cpp
##    src/isometric/unary/compare_scalar_helpers.cpp
##    src/isometric/unary/compare_vector_helpers.cpp
##    src/isometric/unary/boolean_scalar_helpers.cpp
##    src/isometric/unary/boolean_vector_helpers.cpp
    src/isometric/binary/arith_helpers.cpp
#    src/isometric/binary/compare_helpers.cpp
#    src/isometric/binary/boolean_helpers.cpp

    src/sparse/CompressedSparseMatrix.cpp
    src/sparse/FragmentedSparseMatrix.cpp
    src/sparse/secondary_extraction.cpp
    src/sparse/convert_to_compressed_sparse.cpp
    src/sparse/convert_to_fragmented_sparse.cpp
    src/sparse/compress_sparse_triplets.cpp

    src/stats/sums.cpp
    src/stats/variances.cpp
    src/stats/medians.cpp
    src/stats/ranges.cpp
    src/stats/counts.cpp
    src/stats/grouped_medians.cpp
    src/stats/grouped_sums.cpp

    src/utils/parallelize.cpp
    src/utils/wrap_shared_ptr.cpp
    src/utils/SomeNumericArray.cpp
    src/utils/ArrayView.cpp
    src/utils/ConsecutiveOracle.cpp
    src/utils/FixedOracle.cpp
    src/utils/process_consecutive_indices.cpp
)

target_link_libraries(
    libtest
    gtest_main
    gmock_main
    tatami
)

target_compile_definitions(libtest PRIVATE DEBUG=1)
#target_compile_options(libtest PRIVATE -Wall -Wextra -Wpedantic -Werror)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
set(DO_CODE_COVERAGE OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
    set(DO_CODE_COVERAGE ON)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)

# Test custom parallelization during apply.
macro(create_partest target)
    add_executable(
        ${target}
        src/stats/sums.cpp
        src/stats/variances.cpp
        src/stats/medians.cpp
        src/stats/ranges.cpp
        src/stats/counts.cpp
        src/utils/parallelize.cpp
        src/stats/grouped_medians.cpp
        src/stats/grouped_sums.cpp
    )

    target_link_libraries(
        ${target}
        gtest_main
        tatami
    )

    #target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)

    if(DO_COVERAGE)
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif()

    gtest_discover_tests(${target})
endmacro()

create_partest(cuspartest)
target_compile_definitions(cuspartest PRIVATE CUSTOM_PARALLEL_TEST=1)

find_package(OpenMP)
if(OpenMP_FOUND)
    create_partest(omptest)
    target_link_libraries(omptest OpenMP::OpenMP_CXX)
endif()
