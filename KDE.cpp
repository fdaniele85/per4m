//
// Created by Tommaso on 27/03/2024.
//

#include "KDE.h"
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ffp {
    bool icase_compare(const std::string_view s1, const std::string_view s2) {
        auto comparator = [](const char a, const char b) { return std::tolower(a) == std::tolower(b); };

        return (s1.size() == s2.size()) &&
               std::equal(s1.begin(), s1.end(), s2.begin(), comparator);
    }

    double gaussianKernel(double x, double bandwidth) { return exp(-0.5 * pow(x, 2) / pow(bandwidth, 2)) / (sqrt(2 * M_PI) * bandwidth); }

    double epanechnikovKernel(double x, double bandwidth) {
        const auto base = x / bandwidth;
        if (std::abs(x) <= bandwidth) {
            return (3.0 / 4.0) * (1 - pow(base, 2));
        } else {
            return 0;
        }
    }

    double uniformKernel(double x, double bandwidth) { return (0.5 / bandwidth) * (std::abs(x) <= bandwidth); }

    KDE::KDE(const Kernel kernel_type, const double lb) : kernel_type_(kernel_type), lb_(lb) {
        switch (kernel_type_) {
        case Kernel::gaussian: kernel_function_ = &gaussianKernel; break;
        case Kernel::epanechnikov: kernel_function_ = &epanechnikovKernel; break;
        case Kernel::uniform: kernel_function_ = &uniformKernel; break;
        }
    }

    std::vector<double> KDE::kernelDensityEstimate(const std::vector<double> &x_values, const std::vector<double> &data) const {
        const auto n = data.size();
        std::vector<double> result(x_values.size(), 0.0);
        std::vector<double> kernel_values(x_values.size(), 0.0);

        for (uint i = 0; i < n; ++i) {
            for (uint j = 0; j < x_values.size(); ++j) {
                kernel_values[j] = kernel_function_(x_values[j] - data[i], bandwidth_);

                if (kernel_values[j] > 0) {
                    result[j] += kernel_values[j] / (static_cast<double>(n) * bandwidth_);
                }
            }
        }
        return result;
    }

    void KDE::feed_data(const std::vector<double> &data) {
        data_ = data;

        double sigma = 0.0;
        const auto n = static_cast<double>(data_.size());
        const auto accumulation = std::accumulate(data_.begin(), data_.end(), 0.0);
        for (const double val : data_) {
            sigma += pow(val - (accumulation / n), 2);
        }
        sigma = sqrt(sigma / n);
        bandwidth_ = pow((4 * pow(sigma, 5)) / (3 * n), 0.2);
    }

    double KDE::cdf(const double x_value, const int n_queries) const {
        std::vector<double> x_values;

        double current_val = lb_;
        const double delta_val = (x_value - lb_) / n_queries;

        for (int i = 0; i < n_queries; ++i) {
            x_values.push_back(current_val);
            current_val += delta_val;
        }

        // TODO implementare per la massimizzazione

        const std::vector<double> pdf_values = kernelDensityEstimate(x_values, data_);

        double cdf_value = 0.0;
        for (uint i = 1; i < pdf_values.size(); ++i) {
            cdf_value += 0.5 * (pdf_values[i] + pdf_values[i - 1]) * (x_values[i] - x_values[i - 1]);
        }
        return cdf_value;
    }

    Kernel get_kernel(std::string_view kernel) {
        if (icase_compare(kernel, "gaussian")) {
            return Kernel::gaussian;
        }

        if (icase_compare(kernel, "epanechnikov")) {
            return Kernel::epanechnikov;
        }

        if (icase_compare(kernel, "uniform")) {
            return Kernel::uniform;
        }

        throw std::invalid_argument("Kernel not defined");
    }
} // namespace ffp