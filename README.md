# A C++ API for all sorts of matrices 

![Unit tests](https://github.com/LTLA/tatami/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/LTLA/tatami/actions/workflows/doxygenate.yaml/badge.svg)

## Overview

**tatami** is a spiritual successor to the [**beachmat** C++ API](https://github.com/LTLA/beachmat) that provides read access to different matrix representations.
Specifically, applications can use **tatami** to read rows and/or columns of a matrix without any knowledge of the specific matrix representation.
This allows application developers to write a single piece of code that will work seamlessly with different inputs, even if the underlying representation varies at run-time.

**tatami** (and **beachmat** before it) is motivated by analyses of processed genomics data, where matrices are typically interpreted as a collection of row- or column-wise vectors.
Many applications involve looping over rows or columns to compute some statistic or summary - for example, testing for differential expression within each row of the matrix.
**tatami** is largely optimized for this access pattern.

Representations currently supported by **tatami** include:

- Dense row/column major matrices, with user-defined storage modes and containers.
- Compressed sparse row/column matrices, with user-defined storage modes and containers.
- Matrices generated by delayed operations, ~stolen~ inspired by those in [**DelayedArray**](https://github.com/Bioconductor/DelayedArray).

## Quick start

**tatami** is a header-only library so can be easily used by just including the source files and creating the desired matrices:

```cpp
#include "tatami/tatami.h"

std::shared_ptr<tatami::numeric_matrix> mat(new tatami::DenseRowMatrix<double>(nrows, ncols, vals));

// Get the dimensions:
size_t NR = mat->nrow(), NC = mat->ncol();

// Extract a column 'i':
std::vector<double> buffer(NR);
auto ptr = mat->column(i, buffer.data());
ptr[0]; // first element of the column.

// Extract the [5, 10) entries of the 'i'-th column.
auto ptr = mat->column(i, buffer.data(), 5, 10);
```

The key idea here is that, once `mat` is created, the application does not need to worry about the exact format of the matrix referenced by the pointer.
Application developers can then write code that works interchangeably with a variety of different matrix representations.

## Access patterns

The **tatami** API is designed for read-only access - it is not possible to alter the matrix contents via the API.
This is especially relevant for matrices with delayed operations or those referring to remote data stores, where reading the matrix data is trivial but writing is not guaranteed to work.
As a result, supported operations are limited to reading data from the matrix:

- `nrow()` and `ncol()` return the number of rows and columns, respectively.
- `row()` and `column()` return pointers to the start of a (contiguous slice of a) row and column, respectively.
- `sparse_row()` and `sparse_column()` return pointers to the values and indices of the non-zero elements in a row and column, respectively.

```cpp
std::vector<double> ibuffer(NC), vbuffer(NC);
auto indexed = mat->sparse_row(i, vbuffer.data(), ibuffer.data());
for (size_t i = 0; i < indexed.number; ++i) {
    indexed.index[i]; // index of the element
    indexed.value[i]; // value of the element
}
```

For performance-critical sections, it may be desirable to customize the client code based on properties of the matrix.
This is supported with the following methods:

- `sparse()` indicates whether a matrix is sparse.
- `prefer_rows()` indicates whether a matrix is more efficiently access along its rows (e.g., row-major dense matrices).

This allows client developers to design special code paths to take advantage of these properties - 
the [`colsums.cpp`](https://github.com/LTLA/tatami/tree/master/gallery/colsums.cpp) example is particularly demonstrative.

All methods in **tatami** are `const` and thus can be used concurrently.
Any mutable information that needs to persist across API calls is handled by passing a writeable pointer to a `workspace` object to each call.
This can be used to cache information across calls for greater efficiency, e.g., when iterating across consecutive rows or columns of a matrix.

```cpp
auto wrk = mat->new_workspace(true);
std::vector<double> buffer2(NC);
for (size_t i = 0; i < mat->nrow(); ++i) {
    auto ptr = mat->row(i, buffer2.data(), wrk.get());
}
```

Check out the [reference documentation](https://ltla.github.io/tatami) for more details on the available classes and operations.
The [gallery](https://github.com/LTLA/tatami/tree/master/gallery) also contains worked examples for common operations based on row/column traversals.

## API design

Matrix classes named with `snake_case` are virtual and intended for use as interfaces - these cannot be directly constructed.
Matrix classes named with `CamelCase` correspond to actual matrix representations and can be explicitly constructed.
All other functions or non-matrix classes use `snake_case`.

We use templating to define the type of values returned by the interface.
This includes the type of the data (most typically `double`) as well as the type of row/column indices (default `int`, but one could imagine using, e.g., `size_t`).
It is worth noting that the storage type does not need to be the same as the interface type.
For example, developers could store a matrix of small counts as `uint16_t` while returning `double`s for compatibility with downstream mathematical code.

Most of the examples in the **tatami** documentation use a `shared_ptr` to a `tatami::numeric_matrix`, relying on run-time polymorphism to determine the right method to call.
This allows applications to easily handle a variety of different input representations.
Alternatively, applications may use templating to achieve compile-time polymorphism on the different **tatami** subclasses,
but this is rather restrictive without providing obvious performance benefits. 

## Other operations

As previously mentioned, **tatami** is really designed to pull out rows or columns of a matrix, and little else.
Some support is provided for basic statistics in the same vein as the [**matrixStats**](https://github.com/HenrikBengtsson/matrixStats) package:

```cpp
auto colsums = tatami::column_sums(mat);
auto rowvars = tatami::row_variances(mat);
```

**tatami** does not support matrix algebra or decompositions.
For this we typically use [**Eigen**](https://eigen.tuxfamily.org/), effectively trading the diversity of representations for a much more comprehensive suite of operations.
A frequent pattern is to use **tatami** to load the input data, which is usually in a custom format to save memory for large datasets;
process it into a smaller submatrix, e.g., by selecting features of interest in a genome-scale analysis;
and then copy this cheaply into an `Eigen::MatrixXd` for more computationally intensive work.

## Building projects with **tatami** 

If you're using CMake, you just need to add something like this to your `CMakeLists.txt`:

```
include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/LTLA/tatami
  GIT_TAG master # or any version of interest 
)

FetchContent_MakeAvailable(tatami)
```

Then you can link to **tatami** to make the headers available during compilation:

```
# For executables:
target_link_libraries(myexe tatami)

# For libaries
target_link_libraries(mylib INTERFACE tatami)
```

If you're not using CMake, the simple approach is to just copy the files - either directly or with Git submodules - and include their path during compilation with, e.g., GCC's `-I`.

## TODO

- Add bindings for TileDB.
