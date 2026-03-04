#pragma once
#include <vector>
#include <optional>

namespace per4m {

    struct Bandwidth {
        static double scott_1d(const std::vector<double>& x);
        static double silverman_1d(const std::vector<double>& x);

#ifdef PER4M_USE_FFTW
        static double isj_1d(const std::vector<double>& x,
                             const std::optional<std::vector<double>>& weights = std::nullopt);
#endif
    };

} // namespace per4m