# `tatami` unit tests

We use the [GoogleTest](https://google.github.io/googletest/) framework to perform unit testing on the **tatami** libary code.
To run the tests, compile the test code and run the binary with the following commands:

```sh
cd ..
cmake -S . -B build
build/tests/libtest
```

The philosophy here is to test for correct access across multiple scenarios.
To demonstrate, consider the tests for column access from a particular matrix representation.
We want to consider all combinations of the following options:

- Access with and without a supplied workspace.
- Access to the full column or to a slice.
  The slice can be constant across calls or can vary across calls;
  the latter checks that the workspace "memory" does not depend on identical calls.
- Access to consecutive columns, or to every _n_-th column.
  This checks that the workspace "memory" does not rely on consecutive access.
- Access in dense and sparse modes.

Where possible, we use parametrized fixtures to avoid code duplication.
The equivalent test requirements apply to row access.

For more complicated representations, further parameter combinations may be necessary.
This is particularly true of matrices that contain other matrices, e.g., as part of delayed operations.
In such cases, it is important to check that the enclosing matrix supports different backends;
we usually test both dense and sparse backends with row- and column-major storage.

The statistic-computing utilities (column sums, row sums, etc.) are more difficult to test without just duplicating the code therein.
As such, we just verify that each utility gives the same result regardless of the backend (i.e., dense/sparse, row- or column-major).
