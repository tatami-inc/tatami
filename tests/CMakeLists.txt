include(FetchContent)
FetchContent_Declare(
  tatami_test
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_test
  GIT_TAG master
)
FetchContent_MakeAvailable(tatami_test)

include(GoogleTest)
set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
set(DO_CODE_COVERAGE OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(DO_CODE_COVERAGE ON)
endif()

macro(decorate_executable target)
    target_link_libraries(${target} tatami tatami_test)

    target_compile_definitions(${target} PRIVATE DEBUG=1)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)

    if(DO_CODE_COVERAGE)
        target_compile_options(${target} PRIVATE -O0 -g --coverage)
        target_link_options(${target} PRIVATE --coverage)
    endif()

    gtest_discover_tests(${target} DISCOVERY_MODE PRE_TEST)
endmacro()

# Main test executable.
add_executable(
    dense_test
    src/dense/DenseMatrix.cpp
    src/dense/convert_to_dense.cpp
    src/dense/transpose.cpp
)
decorate_executable(dense_test)

add_executable(
    other_test
    src/other/DelayedBind.cpp
    src/other/DelayedTranspose.cpp
    src/other/DelayedCast.cpp
    src/other/ConstantMatrix.cpp
)
decorate_executable(other_test)

add_executable(
    subset_test
    src/subset/DelayedSubset.cpp
    src/subset/DelayedSubsetBlock.cpp
)
decorate_executable(subset_test)

add_executable(
    isometric_unary_test
    src/isometric/unary/DelayedUnaryIsometricOperation.cpp
    src/isometric/unary/arithmetic_vector_helpers.cpp
    src/isometric/unary/arithmetic_scalar_helpers.cpp
    src/isometric/unary/math_helpers.cpp
    src/isometric/unary/compare_scalar_helpers.cpp
    src/isometric/unary/compare_vector_helpers.cpp
    src/isometric/unary/boolean_scalar_helpers.cpp
    src/isometric/unary/boolean_vector_helpers.cpp
    src/isometric/unary/substitute_scalar_helpers.cpp
    src/isometric/unary/substitute_vector_helpers.cpp
)
decorate_executable(isometric_unary_test)

add_executable(
    isometric_other_test
    src/isometric/arithmetic_utils.cpp
)
decorate_executable(isometric_other_test)

add_executable(
    isometric_binary_test
    src/isometric/binary/DelayedBinaryIsometricOperation.cpp
    src/isometric/binary/arithmetic_helpers.cpp
    src/isometric/binary/compare_helpers.cpp
    src/isometric/binary/boolean_helpers.cpp
)
decorate_executable(isometric_binary_test)

add_executable(
    sparse_test
    src/sparse/CompressedSparseMatrix.cpp
    src/sparse/FragmentedSparseMatrix.cpp
    src/sparse/secondary_extraction.cpp
    src/sparse/convert_to_compressed_sparse.cpp
    src/sparse/convert_to_fragmented_sparse.cpp
    src/sparse/compress_sparse_triplets.cpp
)
decorate_executable(sparse_test)

add_executable(
    utils_test
    src/utils/parallelize.cpp
    src/utils/wrap_shared_ptr.cpp
    src/utils/SomeNumericArray.cpp
    src/utils/ArrayView.cpp
    src/utils/ConsecutiveOracle.cpp
    src/utils/FixedOracle.cpp
    src/utils/process_consecutive_indices.cpp
    src/utils/miscellaneous.cpp
    src/utils/Index_to_container.cpp
)
decorate_executable(utils_test)

add_executable(
    cuspar_test
    src/utils/parallelize.cpp
    src/dense/convert_to_dense.cpp
    src/sparse/convert_to_compressed_sparse.cpp
    src/sparse/convert_to_fragmented_sparse.cpp
)
decorate_executable(cuspar_test)
target_compile_definitions(cuspar_test PRIVATE CUSTOM_PARALLEL_TEST=1)
