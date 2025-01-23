//
// Created by daniele on 23/01/25.
//

#pragma once

#include <algorithm>
#include <string_view>
#include <vector>

namespace ffp {
    /**
     * @enum Kernel
     * @brief Enumeration of kernel types for KDE.
     */
    enum class Kernel {
        gaussian,     ///< Gaussian kernel
        epanechnikov, ///< Epanechnikov kernel
        uniform       ///< Uniform kernel
    };

    /**
     * @brief Retrieves the kernel type from a string representation.
     * @param kernel A string view representing the kernel type.
     * @return The corresponding Kernel enumeration value.
     */
    Kernel get_kernel(std::string_view kernel);

    /**
     * @class KDE
     * @brief Class for performing Kernel Density Estimation (KDE).
     */
    class KDE {
    public:
        /**
         * @brief Constructor for KDE.
         * @param kernel_type The type of kernel to use.
         * @param lb The lower bound for the data (default is 0).
         */
        explicit KDE(Kernel kernel_type, double lb = 0);

        /**
         * @brief Feeds data to the KDE.
         * @param data A vector of data points.
         */
        void feed_data(const std::vector<double> &data);

        /**
         * @brief Estimates the kernel density for given x values.
         * @param x_values A vector of x values where the density is estimated.
         * @param data A vector of data points.
         * @return A vector of estimated density.
         */
        [[nodiscard]] std::vector<double> kernelDensityEstimate(const std::vector<double> &x_values, const std::vector<double> &data) const;

        /**
         * @brief Computes the cumulative distribution function (CDF) for a given x value.
         * @param x_value The x value for which the CDF is computed.
         * @param n_queries The number of queries to use in the computation.
         * @return The CDF value.
         */
        [[nodiscard]] double cdf(double x_value, int n_queries) const;

    private:
        double bandwidth_{0.5};                    ///< Bandwidth for the kernel.
        Kernel kernel_type_{Kernel::epanechnikov}; ///< Type of kernel used.

        // TODO: Usare std::queue per mantenere solo gli ultimi n_sample valori e stimare sempre su quelli?
        std::vector<double> data_; ///< Data points for KDE.
        double lb_;                ///< Lower bound for the data.

        std::function<double(double, double)> kernel_function_; ///< Kernel function.
    };
} // namespace ffp
