include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    libtest 
    src/base/DenseMatrix.cpp
    src/base/CompressedSparseMatrix.cpp
    src/base/DelayedBind.cpp
    src/base/DelayedSubset.cpp
    src/base/DelayedSubsetBlock.cpp
    src/base/DelayedTranspose.cpp
    src/base/arith_vector_helpers.cpp
    src/base/arith_scalar_helpers.cpp
    src/base/math_helpers.cpp
    src/base/compress_sparse_triplets.cpp
    src/stats/sums.cpp
    src/stats/variances.cpp
    src/stats/medians.cpp
    src/data/data_sparse.cpp
    src/data/data_triangular.cpp
)

target_link_libraries(
    libtest
    gtest_main
    tatami
)

add_executable(
    exttest 
    src/ext/MatrixMarket.cpp
)

target_link_libraries(
    exttest
    gtest_main
    tatami
)

include(GoogleTest)
gtest_discover_tests(libtest)

gtest_discover_tests(exttest)

add_test(NAME tatami_tests COMMAND libtest exttest)
