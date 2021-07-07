# Gallery of `tatami` examples

All of the examples can be compiled with:

```sh
cd ..
cmake -DBUILD_TESTING=OFF -DBUILD_TATAMI_GALLERY=ON -S . -B build
cmake --build build
```

The resulting binaries will be available in the `build/gallery` subdirectory,
named after the C++ source file from which they were generated.

- `colsums.cpp`: use column (and eventually also row) extraction methods for a typical `tatami::numeric_matrix` to compute the column sums.
- `parallel.cpp`: compute row sums in parallel with OpenMP.
- `char2double.cpp`: store integers as `char` to save memory, but return them as `double`s for compatibility with downstream code.
- `sparse_workspace.cpp`: compare sparse matrix access speeds with and without a workspace.
