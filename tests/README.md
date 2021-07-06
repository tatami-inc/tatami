# **tatami** unit tests

We use the [GoogleTest](https://google.github.io/googletest/) framework to perform unit testing on the **tatami** libary code.
To run the tests, compile the test code and run the binary:

```sh
cd ..
cmake -S . -B build
build/tests/libtest
```

For matrix representations, the philosophy is to check access under multiple scenarios.
For column access, we want to consider all combinations of:

- Access with and without a supplied workspace.
- Access to the full column or to a slice.
  The slice can be constant across calls or can vary across calls;
  the latter checks that the workspace "memory" does not depend on identical calls.
- Access to consecutive columns, or to every _n_-th column.
  This checks that the workspace "memory" does not rely on consecutive access.
- Access in dense and sparse modes.

The equivalent requirements apply to row access.
Where possible, we use parametrized fixtures to avoid code duplication.

For more complicated representations, further parameter combinations may be necessary.
This is particularly true of matrices that contain other matrices, e.g., as part of delayed operations.
In such cases, it is important to check that the enclosing matrix supports different backends;
we usually test both dense and sparse backends with row- and column-major storage.

The statistic-computing utilities (column sums, row sums, etc.) are more difficult to test without just duplicating the code therein.
As such, we just verify that each utility gives the same result regardless of the backend (i.e., dense/sparse, row- or column-major).
