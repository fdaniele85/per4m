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

## How to integrate probabilistic stopping in an existing algorithm

`kde_stop::ProbabilisticStop` is designed to be added to an existing optimization loop with minimal changes:

1. create one stopper before the main loop
2. feed it objective values with `add(cost)`
3. check `stop()` after each meaningful observation
4. keep your usual deterministic limits, such as time limit or maximum iterations, as hard safeguards

The stopping rule is intended for minimization. It estimates the probability of observing a solution that improves the best value seen so far by at least `improve_pct`. When that estimated probability falls below `threshold`, `stop()` returns `true`.

```cpp
#include <kde_stop/ProbabilisticStop.h>

kde_stop::ProbabilisticStop stopper(
    0.01,                              // threshold: stop below 1% improvement probability
    0.001,                             // improve_pct: require at least 0.1% relative improvement
    kde_stop::Kernel::gaussian,         // KDE kernel
    known_objective_lower_bound,        // lower bound on the observed cost scale
    50,                                // number of recent observations used by the KDE
    1024,                              // number of CDF integration queries
    kde_stop::BandwidthType::silverman  // bandwidth rule
);
```

Parameter notes:

- `threshold`: smaller values make the criterion more conservative
- `improve_pct`: the relative improvement worth waiting for, for example `0.001` for 0.1%
- `known_objective_lower_bound`: a lower bound for the cost values passed to `add`, and below the target values being integrated by the KDE
- the fifth constructor argument is the sliding window size used by the KDE
- `number_of_queries` controls the numerical integration resolution; values around `512` to `2048` are a practical starting point
- `stop()` remains `false` until enough observations have been added to fill the first KDE window
- the internal improvement target is `best - improve_pct * best`; for the usual relative-improvement interpretation, feed positive minimization costs
- for maximization, pass a transformed minimization value, for example `-score`, and choose a consistent lower bound for that transformed scale

### Example: GRASP

In a GRASP, each iteration typically builds a randomized solution and improves it with local search. A natural integration point is after local search, where each iteration contributes one locally optimal cost sample.

```cpp
#include <utility>
#include <kde_stop/ProbabilisticStop.h>

Solution run_grasp(const Instance& instance, int max_iterations, double lower_bound) {
    constexpr double stop_probability = 0.01;
    constexpr double required_improvement = 0.001;
    constexpr int window_size = 50;
    constexpr int kde_queries = 1024;

    kde_stop::ProbabilisticStop stopper(
        stop_probability,
        required_improvement,
        kde_stop::Kernel::gaussian,
        lower_bound,
        window_size,
        kde_queries,
        kde_stop::BandwidthType::silverman
    );

    Solution best;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        Solution candidate = construct_randomized_solution(instance);
        candidate = local_search(instance, std::move(candidate));

        const double cost = candidate.cost();
        if (!best.is_valid() || cost < best.cost()) {
            best = candidate;
        }

        stopper.add(cost);
        if (stopper.stop()) {
            break;
        }
    }

    return best;
}
```

This feeds the KDE with the distribution of GRASP iteration outcomes, not only with the incumbent trajectory. That is usually preferable because the best-so-far value is monotone and can quickly become a long plateau.

### Example: ALNS

In an ALNS, the stopper can be checked after each iteration, exactly as in the GRASP example. Feed the objective value produced by the current ALNS step, then stop when the estimated probability of obtaining the required improvement becomes too small.

```cpp
#include <utility>
#include <kde_stop/ProbabilisticStop.h>

Solution run_alns(const Instance& instance, int max_iterations, double lower_bound) {
    constexpr double stop_probability = 0.005;
    constexpr double required_improvement = 0.0005;
    constexpr int window_size = 30;
    constexpr int kde_queries = 1024;

    kde_stop::ProbabilisticStop stopper(
        stop_probability,
        required_improvement,
        kde_stop::Kernel::epanechnikov,
        lower_bound,
        window_size,
        kde_queries,
        kde_stop::BandwidthType::silverman
    );

    Solution current = initial_solution(instance);
    Solution best = current;

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        auto destroy = select_destroy_operator();
        auto repair = select_repair_operator();

        Solution candidate = repair(instance, destroy(instance, current));

        if (accept(candidate, current, iteration)) {
            current = std::move(candidate);
        }

        if (current.cost() < best.cost()) {
            best = current;
        }

        update_operator_scores(current, best);
        update_operator_weights();

        stopper.add(current.cost());
        if (stopper.stop()) {
            break;
        }
    }

    return best;
}
```

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
