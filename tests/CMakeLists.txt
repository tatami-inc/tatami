include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/afd902e992b720d1b3e106bc5e425a5768872265.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Main test executable.
add_executable(
    libtest 

    src/dense/DenseMatrix.cpp

    src/other/DelayedBind.cpp
    src/other/DelayedTranspose.cpp
    src/other/DelayedCast.cpp

    src/subset/DelayedSubset.cpp
    src/subset/DelayedSubsetBlock.cpp

    src/isometric/unary/arith_vector_helpers.cpp
    src/isometric/unary/arith_scalar_helpers.cpp
    src/isometric/unary/math_helpers.cpp
    src/isometric/unary/compare_scalar_helpers.cpp
    src/isometric/unary/compare_vector_helpers.cpp
    src/isometric/unary/boolean_scalar_helpers.cpp
    src/isometric/unary/boolean_vector_helpers.cpp
    src/isometric/binary/arith_helpers.cpp
    src/isometric/binary/compare_helpers.cpp
    src/isometric/binary/boolean_helpers.cpp

    src/sparse/CompressedSparseMatrix.cpp
    src/sparse/FragmentedSparseMatrix.cpp
    src/sparse/SemiCompressedSparseMatrix.cpp
    src/sparse/SparseSecondaryExtractorCore.cpp

    src/stats/sums.cpp
    src/stats/variances.cpp
    src/stats/medians.cpp
    src/stats/ranges.cpp
    src/stats/parallelize.cpp

    src/utils/wrap_shared_ptr.cpp
    src/utils/SomeNumericArray.cpp
    src/utils/ArrayView.cpp
    src/utils/convert_to_dense.cpp
    src/utils/convert_to_sparse.cpp
    src/utils/bind_intersection.cpp
    src/utils/compress_sparse_triplets.cpp
    src/utils/Oracles.cpp
    src/utils/process_consecutive_indices.cpp
)

target_link_libraries(
    libtest
    gtest_main
    gmock_main
    tatami
)

# Test custom parallelization during apply.
add_executable(
    cuspartest 
    src/stats/sums.cpp
    src/stats/variances.cpp
    src/stats/medians.cpp
    src/stats/ranges.cpp
    src/stats/parallelize.cpp
)

target_link_libraries(
    cuspartest 
    gtest_main
    tatami
)

# Checking whether to add OpenMP support. This is turned off
# by default to make it easier to debug test failures.
set(USE_OPENMP OFF CACHE BOOL "Compile with OpenMP support")
if (USE_OPENMP)
    find_package(OpenMP)
    target_link_libraries(libtest OpenMP::OpenMP_CXX)
    target_link_libraries(cuspartest OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(libtest PRIVATE DEBUG=1)
target_compile_definitions(cuspartest PRIVATE CUSTOM_PARALLEL_TEST=1)

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
endif()

# Making the tests discoverable.
include(GoogleTest)
gtest_discover_tests(libtest)
gtest_discover_tests(cuspartest)
