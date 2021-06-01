# Gallery of `tatami` examples

## `colsums.cpp`

Use column (and eventually also row) extraction methods for a typical `tatami::numeric_matrix` to compute the column sums.

```sh
g++ -std=c++1z -I../include colsums.cpp
```

## `parallel.cpp`

Compute row sums in parallel with OpenMP.

```sh
g++ -std=c++1z -I../include -fopenmp parallel.cpp
```

## `char2double.cpp`

Store integers as `char` to save memory, but return them as `double`s for compatibility with downstream code.

```sh
g++ -std=c++1z -I../include char2double.cpp
```
