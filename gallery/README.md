# Gallery of `tatami` examples

All of the examples can be compiled with CMake using the commands below:

```sh
cd ..
cmake -DBUILD_TESTING=OFF -DBUILD_TATAMI_GALLERY=ON -S . -B build
cmake --build build
```

This generates executables in the `build/gallery` subdirectory:

- `colsums`: use various extraction methods for a typical `tatami::numeric_matrix` to compute the column sums.
- `parallel`: compute row sums in parallel with OpenMP.
- `char2double`: store integers as `char` to save memory, but return them as `double`s for downstream use.
- `sparse_workspace`: compare sparse matrix access speeds with and without a workspace.

Each executable is named after the C++ source file from which they were generated.
Each file contains some commentary explaining the rationale behind each example.

