# A C++ API for all sorts of matrices 

![Unit tests](https://github.com/tatami-inc/tatami/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/tatami-inc/tatami/actions/workflows/doxygenate.yaml/badge.svg)
[![Codecov](https://codecov.io/gh/tatami-inc/tatami/branch/master/graph/badge.svg?token=Z189ORCLLR)](https://codecov.io/gh/tatami-inc/tatami)

## Overview

**tatami** is a spiritual successor to the [**beachmat** C++ API](https://github.com/tatami-inc/beachmat) that provides read access to different matrix representations.
Specifically, applications can use **tatami** to read rows and/or columns of a matrix without any knowledge of the specific matrix representation.
This allows application developers to write a single piece of code that will work seamlessly with different inputs, even if the underlying representation varies at run-time.
In particular, **tatami** is motivated by analyses of processed genomics data, where matrices are often interpreted as a collection of row- or column-wise vectors.
Many applications involve looping over rows or columns to compute some statistic or summary - for example, testing for differential expression within each row of the matrix.
**tatami** aims to optimize this access pattern across a variety of different matrix representations, depending on how the data is provided to the application.

## Quick start

**tatami** is a header-only library, so it can be easily used by just `#include`ing the relevant source files:

```cpp
#include "tatami/tatami.hpp"

int nrows = 10;
int ncols = 20;
std::vector<double> vals(nrows * ncols);

// 'double' is the data type, 'int' is the row/column index type.
std::shared_ptr<tatami::Matrix<double, int> > mat(
    new tatami::DenseMatrix<double, int, std::vector<double> >(
        nrows,
        ncols,
        vals,
        /* row_major = */ true
    )
);

// Get the dimensions:
int NR = mat->nrow(), NC = mat->ncol();

// Extract the 'i'-th column.
auto extractor = tatami::make_DenseWrapper(*mat, /* row = */ false, {});
auto ptr = extractor.fetch(i);
ptr[0]; // 1st element of the 'i'-th column.

// Extract a block of rows of the 'i'-th column.
auto sliced_extractor = tatami::make_DenseWrapper(*mat, /* row = */ false, 5, 7, {})
auto sliced_ptr = sliced_extractor->fetch(i);
sliced_ptr[0]; // 6th element of the 'i'-th column.

// Extract the first few odd elements of the 'i'-th column.
tatami::VectorPtr<int> indices(new std::vector<int>{ 1, 3, 5, 7, 9 });
auto indexed_extractor = tatami::make_DenseWrapper(*mat, /* row = */ false, indices, {})
auto indexed_ptr = indexed_extractor->fetch(i);
indixed_ptr[0]; // 2nd element of the 'i'-th column.
```

The key idea here is that, once `mat` is created, the application does not need to worry about the exact format of the matrix referenced by the pointer.
Application developers can write code that works interchangeably with a variety of different matrix representations.

## Creating a `tatami::Matrix`

Users can create an instance of a concrete `tatami::Matrix` subclass by using one of the constructors or the equivalent `make_*` utility:

| Description                                | Class                             |
|--------------------------------------------|-----------------------------------|
| Dense matrix                               | `DenseMatrix`                     |
| Compressed sparse matrix                   | `CompressedSparseMatrix`          |
| List of lists sparse matrix                | `FragmentedSparseMatrix`          |
| Delayed isometric unary operation          | `DelayedUnaryIsometricOperation`  |
| Delayed isometric binary operation         | `DelayedBinaryIsometricOperation` |
| Delayed combination                        | `DelayedBind`                     |
| Delayed subset                             | `DelayedSubset`                   |
| Delayed transpose                          | `DelayedTranspose`                |
| Delayed cast                               | `DelayedCast`                     |
| Delayed cast                               | `DelayedCast`                     |
| Constant matrix                            | `ConstantMatrix`                  |

For example, to create a compressed sparse matrix from sparse triplet data, we could do:

```cpp
#include "tatami/tatami.hpp"

int NROW = 10, NCOL = 20;
std::vector<double> x; // vector of non-zero values.
std::vector<int> i; // row indices of length equal to 'x'.
std::vector<int> j; // column indices of length equal to 'x'.

auto indptrs = tatami::compress_sparse_triplets(NROW, NCOL, x, i, j, /* csr = */ false);

std::shared_ptr<tatami::Matrix<double, int> > smat(
    new tatami::CompressedSparseMatrix<
        double, // data type
        int, // index type
        std::vector<double>, // type of 'x'
        std::vector<int>, // type of the row indices 'i'
        std::vector<size_t> // type of the pointers 'indptrs'
    >(
        NR, 
        NC, 
        std::move(x), 
        std::move(i), 
        std::move(indptrs),
        /* csr = */ false
    )
);
```

We typically create a `shared_ptr` to a `tatami::Matrix` to leverage run-time polymorphism.
This enables downstream applications to accept many different matrix representations by compiling against the `tatami::Matrix` interface.
Alternatively, applications may use templating to achieve compile-time polymorphism on the different **tatami** subclasses,
but this is rather restrictive without providing obvious performance benefits. 

We use templating to define the type of values returned by the interface.
This includes the type of the data (most typically `double`) as well as the type of row/column indices (default `int`, but one could imagine using, e.g., `size_t`).
It is worth noting that the storage type does not need to be the same as the interface type.
For example, developers could store a matrix of small counts as `uint16_t` while returning `double`s for compatibility with downstream mathematical code.

The delayed operations are ~stolen from~ inspired by those in the [**DelayedArray**](https://github.com/Bioconductor/DelayedArray) package.
Isometric operations are particularly useful as they accommodate matrix-scalar/vector arithmetic and various mathematical operations.
For example, we could apply a sparsity-breaking delayed operation to our sparse matrix `smat` without actually creating a dense matrix:

```cpp
std::shared_ptr<tatami::Matrix<double, int> > mat2(
    new tatami::DelayedUnaryIsometricOperation<
        double, // type of the result of the operation 
        double, // type of the original matrix
        int // row/column index type
    >(
        smat, 
        std::make_shared<
            tatami::DelayedUnaryIsometricAddScalarHelper<
                double, // type of the result of the operation
                double, // type of the original matrix
                int, // row/column index type
                double // type of the scalar
            >
        >(2.0)
    )
);
```

Some libraries in the [**@tatami-inc**](https://github.com/tatami-inc) organization implement further extensions of **tatami**'s interface, 
e.g., for [HDF5-backed matrices](https://github.com/tatami-inc/tatami_hdf5) and [R-based matrices](https://github.com/tatami-inc/tatami_r).

## Extracting matrix contents

Given an abstract `tatami::Matrix`, we create an `Wrapper` instance to actually extract the matrix data.
Each `Wrapper` object can store intermediate data for re-use during iteration through the matrix, which is helpful for some matrix implementations that do not easily support random access.
The `MyopicDenseWrapper::fetch()` method returns a pointer to the row's contents, valid until the next `fetch()` call on the same instance.
For example, to extract dense rows from our `mat`:

```cpp
const int NR = mat->nrow();
const int NC = mat->ncol();

auto wrp = tatami::make_DenseWrapper(*mat, true, {});
for (int r = 0; r < NR; ++r) {
    auto current = wrp.fetch(r);
    // Do something with the 'current' pointer.
}
```

Alternatively, we could extract sparse columns via `MyopicSparseWrapper::fetch()`.
This returns a `tatami::SparseRange` containing pointers to arrays of (structurally non-zero) values and their row indices,
which provides some opportunities for optimization in algorithms that only need to operate on non-zero values.
For example, to extract sparse columns from our `mat`:

```cpp
auto swrp = tatami::make_SparseWrapper(*mat, /* row = */ false, {});
for (int c = 0; c < NC; ++c) {
    auto current = swrp.fetch(c);
    // 'current' is valid until the next 'swrp.fetch()' call.
    current.number; // Number of structural non-zeros
    current.value; // Pointer to the value array
    current.index; // Pointer to the index array
}
```

In both the dense and sparse cases, we can restrict the values that are extracted by `fetch()`.
This provides some opportunities for optimization by avoiding the unnecessary extraction of uninteresting data.
To do so, we specify the start and length of a contiguous block of interest, or we supply a vector containing the indices of the elements of interest:

```cpp
// Get rows [5, 17) from each column.
auto bwrp = tatami::make_SparseWrapper(*mat, false, 5, 12, {}); 

// Get the first few odd-indexed columns from each row.
auto iwrp = tatami::make_DenseWrapper(
    *mat,
    true,
    tatami::VectorPtr<int>(new std::vector<int>{ 1, 3, 5, 7 }),
    {}
);
```

For sparse extraction, an application might not need the values and/or indices of the structural non-zero elements.
We can disable their extraction for greater efficiency:

```cpp
tatami::Options opt;
opt.sparse_extract_value = false;
opt.sparse_extract_index = false;
auto swrp = tatami::make_SparseWrapper(*mat, false, opt);
auto res = set.fetch(0);
res.number; // only the number of structural non-zeros is returned.
res.value == NULL;
res.index == NULL;
```

## Handling different access patterns

In performance-critical sections, it may be desirable to customize the extraction based on properties of the matrix.
This is supported with the following methods:

- `tatami::Matrix::is_sparse()` indicates whether a matrix is sparse.
- `tatami::Matrix::prefer_rows()` indicates whether a matrix is more efficiently accessed along its rows (e.g., row-major dense matrices).

Users can then write dedicated code paths to take advantage of these properties.
For example, we might use different algorithms for dense data, where we don't have to look up indices; and for sparse data, if we can avoid the uninteresting zero values.
Similarly, if we want to compute a row-wise statistic, but the matrix is more efficiently accessed by column according to `prefer_rows()`,
we could iterate on the columns and attempt to compute the statistic in a "running" manner 
(see [`colsums.cpp`](https://github.com/tatami-inc/gallery/tree/master/src/colsums.cpp) for an example).
In the most complex cases, this leads to code like:

```cpp
if (mat->is_sparse()) {
    if (mat->prefer_rows()) {
        auto swrp = tatami::make_SparseWrapper(*mat, true, {});
        // Do compute along sparse rows.
    } else {
        auto swrp = tatami::make_SparseWrapper(*mat, false, {});
        // Do compute along sparse columns.
    }
} else {
    if (mat->prefer_rows()) {
        auto dwrp = tatami::make_DenseWrapper(*mat, true, {});
        // Do compute along dense rows.
    } else {
        auto dwrp = tatami::make_DenseWrapper(*mat, false, {});
        // Do compute along dense columns.
    }
}
```

Of course, this assumes that it is possible to provide sparse-specific optimizations as well as running calculations for the statistic of interest.
In most cases, only a subset of the extraction patterns are actually feasible so special code paths would not be beneficial.

## Supporting parallelization

The mutable nature of an `Wrapper` instance means that the `fetch()` calls themselves are not `const`.
This means that the same instance cannot be safely re-used across different threads as each call to `fetch()` will modify the instance's contents.
Fortunately, the solution is simple - just create a separate `Wrapper` (and the associated buffers) for each thread.
With OpenMP, this looks like:

```cpp
#pragma omp parallel num_threads(nthreads);
{
    auto wrk = tatami::make_DenseWrapper(*mat, true, {});

    #pragma omp for
    for (int r = 0; r < NR; ++r) {
        auto ptr = wrk->fetch(r);
        // Do something in each thread.
    }
}
```

Users may also consider using the `tatami::parallelize()` function, which accepts a function with the range of jobs (in this case, rows) to be processed in each thread.
This automatically falls back to the standard `<thread>` library if OpenMP is not available.
Applications can also set the `TATAMI_CUSTOM_PARALLEL` macro to override the parallelization scheme in all `tatami::parallelize()` calls.

```cpp
tatami::parallelize([&](int thread, int start, int length) -> void {
    auto wrk = tatami::make_DenseWrapper(*mat, true, {});
    for (int r = 0; r < length; ++r) {
        auto ptr = wrk.fetch(r + start);
        // Do something in each thread.
    }
}, NR, nthreads);
```

## Defining an oracle

When constructing a `Wrapper`, users can supply an `Oracle` that specifies the sequence of rows/columns to be accessed.
Knowledge of the future access pattern enables optimizations in some `Matrix` implementations, 
e.g., file-backed matrices can reduce the number of disk reads by pre-fetching the right data for future accesses.
The most obvious use case involves accessing consecutive rows/columns:

```cpp
auto wrk = tatami::make_DenseWrapper(*mat, true, tatami::make_ConsecutiveOracle(0, NR), {});
for (int r = 0; r < NR; ++r) {
    // No need to specify the index, as 'wrk' already knows the
    // sequence of indices as prediced by the oracle.
    auto ptr = wrk.next();
}
```

Alternatively, we can use the `FixedVectorOracle` class with an array of row/column indices to be extracted:

```cpp
// Extracting columns in an arbitrary predefined order:
std::vector<int> predictions{ 9, 1, 2, 3, 13, 9, 1 };
auto foracle = std::make_shared<tatami::FixedVectorOracle<int> >(std::move(predictions));
auto swrk = tatami::make_DenseWrapper(*mat, false, std::move(foracle), {});
for (std::size_t p = 0; p < predictions.size(); ++p) {
    // No need to specify the index to 'next()', as 'swrk' already knows the
    // sequence of indices as prediced by the oracle.
    auto current = swrk.next();
}
```

Advanced users can also define their own `Oracle` subclasses to generate predictions on the fly.

## Using custom buffers

The `Wrapper` class is just a user-friendly wrapper around an virtual `Extractor` class.
The former allocates the temporary buffers that might be used in each `fetch()` operation,
while subclasses of the latter are responsible for the heavy lifting of pulling data out of each matrix representation.
For users who want more control over memory management, the latter can be used directly, e.g., to extract matrix contents directly into a user-specified array:

```cpp
int NR = mat->nrow();
int NC = mat->ncol();

// Equivalent to 'make_DenseWrapper(*mat, true, {})'
auto ext = mat->dense(true, {});
std::vector<double> buffer(NC);
for (int r = 0; r < NR; ++r) {
    auto current = ext->fetch(r, buffer.data());

    // Note that the contents of 'buffer' is not guaranteed to be the same as 'current'.
    // If we want them to be the same, we need to do some extra work:
    tatami::copy_n(current, NC, buffer.data());

    // Now we can do something with 'buffer', which now contains the r-th row.
}

// Equivalent to 'make_SparseWrapper(*mat, false, 8, 20, {})'
auto sext = mat->sparse(false, 8, 20, {});
std::vector<double> vbuffer(20);
std::vector<int> ibuffer(20);
for (int c = 0; c < NC; ++r) {
    auto current = sext->fetch(c, vbuffer.data(), ibuffer.data());

    // Note that contents of 'vbuffer' and 'ibuffer' might not be the same as 'current'.
    // If we want them to be the same, we need to do some extra work:
    tatami::copy_n(current.value, current.number, vbuffer.data());
    tatami::copy_n(current.index, current.number, ibuffer.data());

    // Now we can do something with 'vbuffer' and 'ibuffer',
    // which now contains the structural non-zeros in this block of the c-th column.
}
```

When developing **tatami**-based applications, we suggest compiling the test suite with and without the `TATAMI_DEBUG_FORCE_COPY` macro.
When set, this macro forces **tatami** to copy data into `buffer`, even if it could omit the copy by directly returning pointers to the internal data structures.
This gives our tests an opportunity to check that a valid buffer is actually supplied to `fetch()`.
Conversely, when the macro is unset and a representation like `DenseMatrix` is used in the test suite,
we can check that our application correctly handles scenarios where `fetch()` does not populate the buffers.

## Comments on other operations

As previously mentioned, **tatami** is primarily designed to pull out rows or columns of a matrix.
Some support is provided for computing basic statistics via the [**tatami_stats**](https://github.com/tatami-inc/tatami_stats) library,
in the same vein as the [**matrixStats**](https://github.com/HenrikBengtsson/matrixStats) R package.
Matrix multiplication is similarly implemented via the [**tatami_mult**](https://github.com/tatami-inc/tatami_mult) library.

**tatami** does not currently support more sophisticated matrix operations like decompositions. 
If these are required, we suggest copying data from **tatami** into other frameworks like [**Eigen**](https://eigen.tuxfamily.org/), 
effectively trading the diversity of representations for a more comprehensive suite of operations.
For example, we often use **tatami** to represent the large input datasets in a custom memory-efficient format;
process it into a much smaller submatrix, e.g., by selecting features of interest in a genome-scale analysis;
and then copy the data into a `Eigen::MatrixXd` or `Eigen::SparseMatrix` for the desired operations.
In practice, many standard decompositions do not scale well for large matrices,
so our applications end up using approximate methods like [**ILRBA**](https://github.com/libscran/irlba) that only require a multiplication operator.

<!---
(At this point, it is worth noting that **Eigen** also supports delayed operations via its expression templates.
These are determined at compile-time and are more efficient than **tatami**'s operations.
However, it is also much harder to pass expression templates around an application while preserving the delayed operations.
Each expression is its own type so supporting multiple expressions would require a potentially-combinatorial increase in the number of realized template functions.
Additionally, the expressions do not extend the lifetime of the matrices on which they operate, so if the matrix is destructed before the expression is evaluated, a segfault will occur.
**tatami** avoids these problems and makes it easier to pass around a matrix with delayed operations.)
-->

**tatami** does not directly support modification of the matrix contents.
Instead, "modifications" are performed by adding delayed operations on top of an immutable matrix.
This avoids difficult bugs where the hypothetical modification of matrix contents via one `shared_ptr` affects all references to the same matrix across the application.
Delayed operations are also more appropriate for matrices referring to read-only data sources, e.g., remote stores or files.
That said, if delayed operations are undesirable, we can use functions like `tatami::convert_to_dense()` to realize our modifications into a new matrix instance. 

## Building projects 

### CMake with `FetchContent`

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```cmake
include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(tatami)
```

Then you can link to **tatami** to make the headers available during compilation:

```cmake
# For executables:
target_link_libraries(myexe tatami)

# For libaries
target_link_libraries(mylib INTERFACE tatami)
```

By default, this will use `FetchContent` to fetch all external dependencies.
Applications are advised to pin the versions of these dependencies themselves - see [`extern/CMakeLists.txt`](extern/CMakeLists.txt) for suggested (minimum) versions.
If you want to install them manually, use `-DTATAMI_FETCH_EXTERN=OFF`.

### CMake with `find_package()`

You can install the library by cloning a suitable version of this repository and running the following commands:

```sh
mkdir build && cd build
cmake .. -DTATAMI_TESTS=OFF
cmake --build . --target install
```

Then you can use `find_package()` as usual:

```cmake
find_package(tatami_tatami CONFIG REQUIRED)
target_link_libraries(mylib INTERFACE tatami::tatami)
```

Again, this will use `FetchContent` to fetch all external dependencies, so see advice above.

### Manual

If you're not using CMake, the simple approach is to just copy the files the `include/` subdirectory - 
either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.
This also requires the external dependencies listed in [`extern/CMakeLists.txt`](extern/CMakeLists.txt). 

## Links

Check out the [reference documentation](https://tatami-inc.github.io/tatami) for more details on each function and class.

The [gallery](https://github.com/tatami-inc/gallery) contains worked examples for common operations based on row/column traversals.

The [**tatami_stats**](https://github.com/tatami-inc/tatami_stats) repository computes some common statistics on **tatami** matrices.

The [**tatami_hdf5**](https://github.com/tatami-inc/tatami_hdf5) repository contains **tatami** bindings for HDF5-backed matrices.

The [**tatami_r**](https://github.com/tatami-inc/tatami_r) repository contains **tatami** bindings for matrix-like objects in R.

The [**beachmat**](https://github.com/tatami-inc/beachmat) package vendors the **tatami** headers for easy use by other R packages.
