# A C++ API for all sorts of matrices 

![Unit tests](https://github.com/LTLA/tatami/actions/workflows/run-tests.yaml/badge.svg)
![Documentation](https://github.com/LTLA/tatami/actions/workflows/doxygenate.yaml/badge.svg)

## Overview

**tatami** is a spiritual successor to the [**beachmat** C++ API](https://github.com/LTLA/beachmat) that provides read access to different matrix representations.
Specifically, applications can use **tatami** to read rows and/or columns of a matrix without any knowledge of the specific matrix representation.
This allows application developers to write a single piece of code that will work seamlessly with different inputs, even if the underlying representation varies at run-time.

**tatami** (and **beachmat** before it) is motivated by analyses of processed genomics data, where matrices are typically interpreted as a collection of row- or column-wise vectors.
Many applications involve looping over rows or columns to compute some statistic or summary - for example, testing for differential expression within each row of the matrix.
**tatami** is largely designed around optimizing this access pattern.

Representations currently supported by **tatami** include:

- Dense row/column major matrices, with user-defined storage modes and containers.
- Compressed sparse row/column matrices, with user-defined storage modes and containers.
- Matrices generated by delayed operations, ~stolen~ inspired by those in [**DelayedArray**](https://github.com/Bioconductor/DelayedArray).

## Quick start

**tatami** is a header-only library so can be easily used by just including the source files:

```cpp
#include "tatami/tatami.h"
```

The key data structure is a shared pointer to a `typed_matrix` class.
This is most commonly a `numeric_matrix`, where the values are specified as `double`s and the row/column indices are specified as `int`s.
Here, we fill up a dense row-major matrix:

```cpp
std::vector<double> vals(50);
std::iota(vals.begin(), vals.end(), 0.0);
std::shared_ptr<tatami::numeric_matrix> mat(new tatami::DenseRowMatrix<double>(10, 5, vals));
```

We can then call a variety of methods without worrying about the exact class of the matrix referenced by `mat`.
For example, to get the number of rows and columns:

```cpp
size_t NR = mat->nrow(), NC = mat->ncol();
```

To extract a row or column, we use the `column()` method to pull out a dense vector.
This requires the caller to supply a buffer though it may not actually be used if the values are already contiguous in memory in the underlying representation.

```cpp
std::vector<double> buffer(NR);
const double* ptr = mat->column(i, buffer.data());
```

And that's it. 

## Concepts

We make heavy use of polymorphism to determine the right method to call at run-time.
This motivates the use of a `shared_ptr` to a `tatami::numeric_matrix` rather than a more conventionally constructed object.
In the example above, after the construction of `mat`, downstream code does not have to care about the exact representation of the matrix data;
this allows clients to write one piece of C++ code that works immediately with a variety of input matrices.

The API is designed for read-only access - it is not possible to alter the matrix contents via the API.
(This is especially true for matrices where the data is remotely hosted, e.g., on AWS S3.)
All methods in the API are similiarly `const` and thus are safe for concurrent use.
Any information that needs to persist across API calls is handled by passing a writeable pointer to a `workspace` object to each call.

For performance-critical sections, it may be desirable to customize the client code based on properties of the matrix.
**tatami** offers the `sparse()` and `prefer_rows()` methods that indicate whether a matrix is sparse and if it prefers extraction on the rows (e.g., for row-major matrices).
This allows client developers to design special code paths to take advantage of these properties - the [`colsums.cpp`](gallery/colsums.cpp) example is particularly demonstrative.

We use templating to define the type of values returned by the interface.
This includes the type of the data (most typically `double`) as well as the type of row/column indices (default `int`, but one could imagine using, e.g., `size_t`).
It is worth noting that the storage type does not need to be the same as the interface type.
For example, developers could store a matrix of small counts as `uint16_t` while returning `double`s for compatibility with downstream mathematical code.

Matrix classes named with `snake_case` are virtual and intended for use as interfaces - these cannot be directly constructed.
Matrix classes named with `CamelCase` correspond to actual matrix representations and can be explicitly constructed.
All other functions or non-matrix classes use `snake_case`.

**tatami** does not support more complex operations like decompositions or matrix algebra.
For this we typically use [**Eigen**](https://eigen.tuxfamily.org/), effectively trading the diversity of representations for a much more comprehensive suite of operations.
A frequent pattern is to use **tatami** to load and preprocess the input data into a smaller submatrix that can be cheaply copied into an `Eigen::MatrixXd` for further work.

## Including **tatami** in a project

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

## Documentation

The [gallery](https://github.com/LTLA/tatami/tree/master/gallery) contains worked examples for common operations based on row/column traversals.

The [`include`](https://github.com/LTLA/tatami/tree/master/include) subdirectory contains further instructions on use, particularly for extensions.

The [reference documentation](https://ltla.github.io/tatami) for the API is generated automatically with Doxygen.

## Running unit tests

We use the [GoogleTest](https://github.com/google/googletest) unit testing framework.
Anyone who is curious can build and execute the tests themselves with:

```sh
cmake -S . -B build
cmake --build build
cd build
make tests
```

## TODO

- Add bindings for TileDB.