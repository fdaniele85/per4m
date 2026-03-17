//
// Created by daniele on 23/01/25.
//

#pragma once

#include "detail/CircularBuffer.h"
#include <functional>
#include <string_view>

namespace per4m {
    /// @enum Kernel
    /// @brief Enumeration of kernel types for KDE.
    enum class Kernel {
        gaussian,     ///< Gaussian kernel
        epanechnikov, ///< Epanechnikov kernel
        uniform       ///< Uniform kernel
    };

    enum class BandwidthType {
        silverman, ///< Silverman's rule of thumb for bandwidth selection
#ifdef PER4M_USE_FFTW
        isj ///< Improved Sheather & Jones method for bandwidth selection
#endif
    };

    /// @brief Retrieves the kernel type from a string representation.
    /// @param kernel A string view representing the kernel type.
    /// @return The corresponding Kernel enumeration value.
    Kernel get_kernel(std::string_view kernel);

    /// @param bandwidth_type A string view representing the bandwidth selection method.
    /// @return The corresponding BandwidthType enumeration value.
    BandwidthType get_bandwidth_type(std::string_view bandwidth_type);

    /// @class KDE
    /// @brief Class for performing Kernel Density Estimation (KDE).
    class KDE {
    public:
        /// @brief Constructor for KDE.
        /// @param kernel_type The type of kernel to use.
        /// @param size The size of the data that will be fed to the KDE.
        /// @param n_queries The number of queries to use in the computation.
        /// @param bandwidth_type The method to be used for bandwidth selection.
        /// @param lb The lower bound for the data (default is 0).
        explicit KDE(Kernel kernel_type, int size, int n_queries, BandwidthType bandwidth_type, double lb = 0);

        /// @brief Destructor for KDE.
        ~KDE() = default;

        /// @brief Feeds data to the KDE.
        /// @param data Data points.
        void feed_data(const detail::CircularBuffer &data);

        /// @brief Estimates the cumulative distribution function (CDF) for a given x value.
        /// @param x_value The x value for which the CDF is computed.
        /// @return The CDF value.
        [[nodiscard]] double estimate(double x_value);

    private:
        double bandwidth_{0.5};                                 ///< Bandwidth for the kernel.
        Kernel kernel_type_{Kernel::epanechnikov};              ///< Type of kernel used.
        BandwidthType bandwidth_type_{BandwidthType::silverman}; ///< Type of bandwidth selection method used.

        std::vector<double> data_; ///< Data points for KDE.
        int size_;                 ///< Size of the data.

        std::vector<double> x_values_;      ///< X values for the KDE.
        std::vector<double> pdf_values_;    ///< Result of the KDE.
        std::vector<double> kernel_values_; ///< Kernel values for the KDE.
        int n_queries_;                     ///< Number of queries.

        double lb_; ///< Lower bound for the data.

        std::function<double(double, double)> kernel_function_{}; ///< Kernel function.

        /// @brief Estimates the kernel density for given x values.
        void kernel_density_estimate();
    };
} // namespace per4m
