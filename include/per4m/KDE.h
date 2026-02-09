//
// Created by daniele on 23/01/25.
//

#pragma once

#include "CircularBuffer.h"
#include <functional>
#include <string_view>

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
         * @param size The size of the data that will be fed to the KDE.
         * @param n_queries The number of queries to use in the computation.
         * @param lb The lower bound for the data (default is 0).
         */
        explicit KDE(Kernel kernel_type, int size, int n_queries, double lb = 0);

        /// @brief Destructor for KDE.
        ~KDE();

        /**
         * @brief Feeds data to the KDE.
         * @param data Data points.
         */
        void feed_data(const CircularBuffer &data);

        /// @brief Estimates the kernel density for given x values.
        void kernelDensityEstimate() const;

        /**
         * @brief Estimates the cumulative distribution function (CDF) for a given x value.
         * @param x_value The x value for which the CDF is computed.
         * @return The CDF value.
         */
        [[nodiscard]] double estimate(double x_value) const;

    private:
        double bandwidth_{0.5};                    ///< Bandwidth for the kernel.
        Kernel kernel_type_{Kernel::epanechnikov}; ///< Type of kernel used.

        double *data_; ///< Data points for KDE.
        int size_;     ///< Size of the data.

        double *x_values_; ///< X values for the KDE.
        double *pdf_values_;  ///< Result of the KDE.
        double *kernel_values_; ///< Kernel values for the KDE.
        int n_queries_;    ///< Number of queries.

        double lb_;                ///< Lower bound for the data.

        std::function<double(double, double)> kernel_function_{}; ///< Kernel function.
    };
} // namespace ffp
