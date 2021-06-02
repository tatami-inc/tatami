cmake_minimum_required(VERSION 3.14)
project(my_project)

set(CMAKE_CXX_STANDARD 17)

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
    src/base/DelayedSubsetOp.cpp
    src/base/arith_vector_helpers.cpp
    src/base/arith_scalar_helpers.cpp
    src/base/math_helpers.cpp
    src/base/compress_sparse_triplets.cpp
    src/data_sparse.cpp
)

target_link_libraries(
    libtest
    gtest_main
)

include_directories(../include)

include(GoogleTest)
gtest_discover_tests(libtest)
