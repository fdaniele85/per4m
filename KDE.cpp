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

        return (s1.size() == s2.size()) && std::equal(s1.begin(), s1.end(), s2.begin(), comparator);
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

    KDE::KDE(const Kernel kernel_type, const int size, const int n_queries, const double lb)
        : kernel_type_(kernel_type),
          data_(new double[size]),
          size_(size),
          x_values_(new double[n_queries]),
          pdf_values_(new double[n_queries]),
          kernel_values_(new double[n_queries]),
          n_queries_(n_queries),
          lb_(lb) {
        switch (kernel_type_) {
        case Kernel::gaussian: kernel_function_ = &gaussianKernel; break;
        case Kernel::epanechnikov: kernel_function_ = &epanechnikovKernel; break;
        default: kernel_function_ = &uniformKernel; break;
        }
    }

    void KDE::kernelDensityEstimate() const {
        for (int i = 0; i < n_queries_; ++i) {
            pdf_values_[i] = 0.0;
            kernel_values_[i] = 0.0;
        }

        for (uint i = 0; i < size_; ++i) {
            for (uint j = 0; j < n_queries_; ++j) {
                kernel_values_[j] = kernel_function_(x_values_[j] - data_[i], bandwidth_);

                if (kernel_values_[j] > 0) {
                    pdf_values_[j] += kernel_values_[j] / (static_cast<double>(size_) * bandwidth_);
                }
            }
        }
    }

    void KDE::feed_data(const CircularBuffer &data) {
        double accumulation = 0.0;
        double sigma = 0.0;

        for (int i = 0; i < size_; ++i) {
            data_[i] = data[i];
            accumulation += data[i];
        }

        for (int i = 0; i < size_; ++i) {
            sigma += pow(data_[i] - (accumulation / size_), 2);
        }
        sigma = sqrt(sigma / size_);
        bandwidth_ = pow((4 * pow(sigma, 5)) / (3 * size_), 0.2);
    }

    double KDE::estimate(const double x_value) const {
        double current_val = lb_;
        const double delta_val = (x_value - lb_) / n_queries_;

        for (int i = 0; i < n_queries_; ++i) {
            x_values_[i] = current_val;
            current_val += delta_val;
        }

        // TODO implementare per la massimizzazione

        kernelDensityEstimate();

        double cdf_value = 0.0;
        for (uint i = 1; i < n_queries_; ++i) {
            cdf_value += 0.5 * (pdf_values_[i] + pdf_values_[i - 1]) * (x_values_[i] - x_values_[i - 1]);
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

    KDE::~KDE() {
        delete[] data_;
        delete[] x_values_;
        delete[] pdf_values_;
        delete[] kernel_values_;
    }
} // namespace ffp