# kde_stop

`kde_stop` is a C++20 library that provides probabilistic stopping utilities based on kernel density estimation (KDE).

## Citation

If you use `kde_stop` in academic work, please cite the paper describing KDE-STOP:

```bibtex
@article{FeroneFestaPastore2026,
  title   = {Enhancing optimization algorithms with Kernel Density Estimation: A statistical learning strategy for smarter metaheuristics},
  volume  = {194},
  doi     = {10.1016/j.cor.2026.107539},
  journal = {Computers \& Operations Research},
  author  = {Ferone, Daniele and Festa, Paola and Pastore, Tommaso},
  year    = {2026},
  month   = oct,
  pages   = {107539}
}
```

The article is available at: http://dx.doi.org/10.1016/j.cor.2026.107539

## Public API

The supported public API currently consists of:

- `kde_stop::ProbabilisticStop`
- `kde_stop::KDE`
- `kde_stop::Bandwidth`

## Requirements

- CMake 3.26+
- A C++20 compiler
- `pkg-config`
- FFTW3 development files (`fftw3`)

On Debian/Ubuntu-like systems, the required packages are typically:

```bash
sudo apt install cmake pkg-config libfftw3-dev doxygen
```

## Build

```bash
cmake -S . -B build
cmake --build build
```

To install the library:

```bash
cmake --install build --prefix /tmp/kde_stop-install
```

## Using the library in another CMake project

After installation:

```cmake
find_package(kde_stop REQUIRED)
target_link_libraries(my_target PRIVATE kde_stop::kde_stop)
```

The package configuration resolves the exported FFTW3 dependency through `pkg-config`, so the consumer project must also have `pkg-config` and FFTW3 available.

## Using the library with `FetchContent`

If you prefer to vendor `kde_stop` directly from CMake instead of installing it first, you can use `FetchContent`:

```cmake
include(FetchContent)

set(KDE_STOP_THREAD_SAFE ON CACHE BOOL "Enable thread-safe mode in kde_stop")

FetchContent_Declare(kde_stop
    URL https://github.com/fdaniele85/kde_stop/archive/refs/tags/v1.0.0.zip
)

FetchContent_MakeAvailable(kde_stop)

target_link_libraries(my_target PRIVATE kde_stop::kde_stop)
```

Notes:

- with `FetchContent`, you do **not** call `find_package(kde_stop)`; the target is created directly by the fetched project
- `KDE_STOP_THREAD_SAFE` must be set **before** `FetchContent_MakeAvailable(kde_stop)` so that `kde_stop` is configured with the desired option
- FFTW3 is still required at configure time, so the host system must provide both `pkg-config` and the `fftw3` development package
- if FFTW3 is not found, CMake configuration still fails because the dependency is `REQUIRED`

## `KDE_STOP_THREAD_SAFE`

The project exposes the CMake option:

```bash
-DKDE_STOP_THREAD_SAFE=ON
```

When enabled, the library adds the public compile definition `KDE_STOP_THREAD_SAFE` and the internal synchronization helpers use `std::mutex`/`std::lock_guard`.

When disabled (default), the same locking helpers become no-op lightweight placeholders, so the code builds without mutex-based synchronization overhead.

Example configure command:

```bash
cmake -S . -B build -DKDE_STOP_THREAD_SAFE=ON
```

## FFTW3 dependency

`kde_stop` currently looks for FFTW3 through `pkg-config`:

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
```

If FFTW3 is found:

- the library links against `PkgConfig::FFTW3`
- the public compile definition `KDE_STOP_USE_FFTW` is enabled
- the FFTW-backed bandwidth code is compiled in

If FFTW3 is **not** found:

- CMake configuration **fails immediately** because the dependency is marked as `REQUIRED`
- no build files are generated for `kde_stop`

At runtime, when FFTW support is available and the ISJ bandwidth computation throws an exception, `kde_stop::KDE` falls back to Silverman's rule.

## Generating documentation with Doxygen

A `Doxyfile` is provided in the project root.

Generate the documentation with:

```bash
doxygen Doxyfile
```

The HTML output is generated under:

```text
docs/doxygen/html/
```

The Doxygen configuration excludes:

- `include/kde_stop/detail/`
- `include/kde_stop/ProbabilisticFilter.h`
- `include/kde_stop/Ribeiro.h`

