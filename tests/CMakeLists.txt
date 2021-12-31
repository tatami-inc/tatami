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
    src/stats/ranges.cpp
    src/stats/apply.cpp
    src/utils/wrap_shared_ptr.cpp
    src/utils/NakedArray.cpp
    src/data/data_sparse.cpp
    src/data/data_triangular.cpp
)

target_link_libraries(
    libtest
    gtest_main
    tatami
)

set(USE_OPENMP OFF CACHE BOOL "Compile with OpenMP support")
if (USE_OPENMP)
    find_package(OpenMP)
    target_link_libraries(libtest OpenMP::OpenMP_CXX)
endif()

add_executable(
    exttest 
    src/ext/MatrixMarketBuffer.cpp
    src/ext/MatrixMarketText.cpp
    src/ext/MatrixMarketGzip.cpp
)

target_link_libraries(
    exttest
    gtest_main
    tatami
)

find_package(ZLIB)
target_link_libraries(exttest ZLIB::ZLIB)
target_compile_definitions(exttest PRIVATE TATAMI_USE_ZLIB=1)

include(CheckIncludeFiles)
check_include_files(filesystem HAVE_CXX_FS)
if (NOT HAVE_CXX_FS) 
    target_link_libraries(exttest stdc++fs) 
endif()

set(CODE_COVERAGE OFF CACHE BOOL "Enable coverage testing")
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(libtest PRIVATE -O0 -g --coverage)
    target_compile_options(exttest PRIVATE -O0 -g --coverage)
    target_link_options(libtest PRIVATE --coverage)
    target_link_options(exttest PRIVATE --coverage)
endif()

include(GoogleTest)
gtest_discover_tests(libtest)

gtest_discover_tests(exttest)

add_test(NAME tatami_tests COMMAND libtest exttest)
