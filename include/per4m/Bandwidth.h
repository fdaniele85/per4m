#pragma once
#include <vector>
#include <optional>

namespace per4m {

    /// \brief A utility class for bandwidth selection methods used in Kernel Density Estimation (KDE).
    struct Bandwidth {
        /// \brief Computes the bandwidth using Silverman's rule of thumb for 1D data.
        /// \param x A vector of data points for which the bandwidth is to be computed.
        /// \return The computed bandwidth value based on Silverman's rule.
        static double silverman_1d(const std::vector<double>& x);

#ifdef PER4M_USE_FFTW
        /// \brief Computes the bandwidth using the Improved Sheather & Jones method for 1D data.
        /// \param x A vector of data points for which the bandwidth is to be computed.
        /// \param weights An optional vector of weights corresponding to the data points. If not provided, equal weights are assumed.
        /// \return The computed bandwidth value based on the ISJ method.
        static double isj_1d(const std::vector<double>& x,
                             const std::optional<std::vector<double>>& weights = std::nullopt);
#endif
    };

} // namespace per4m