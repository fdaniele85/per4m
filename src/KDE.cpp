//
// Created by Tommaso on 27/03/2024.
//

#include <cmath>
#include <iostream>
#include <per4m/Bandwidth.h>
#include <per4m/KDE.h>
#include <stdexcept>

namespace per4m {
    bool icase_compare(const std::string_view s1, const std::string_view s2) {
        auto comparator = [](const char a, const char b) { return std::tolower(a) == std::tolower(b); };

        return (s1.size() == s2.size()) && std::equal(s1.begin(), s1.end(), s2.begin(), comparator);
    }

    double gaussian_kernel(const double x, const double bandwidth) { return exp(-0.5 * pow(x, 2) / pow(bandwidth, 2)) / (sqrt(2 * M_PI) * bandwidth); }

    double epanechnikov_kernel(const double x, const double bandwidth) {
        const auto base = x / bandwidth;
        if (std::abs(x) <= bandwidth) {
            return (3.0 / 4.0) * (1 - pow(base, 2));
        }

        return 0;
    }

    double uniformKernel(const double x, const double bandwidth) { return (0.5 / bandwidth) * (std::abs(x) <= bandwidth); }

    KDE::KDE(const Kernel kernel_type, const int size, const int n_queries, const BandwidthType bandwidth_type, const double lb)
        : kernel_type_(kernel_type),
          bandwidth_type_(bandwidth_type),
          data_(size),
          size_(size),
          x_values_(n_queries),
          pdf_values_(n_queries),
          kernel_values_(n_queries),
          n_queries_(n_queries),
          lb_(lb) {
        switch (kernel_type_) {
        case Kernel::gaussian: kernel_function_ = &gaussian_kernel; break;
        case Kernel::epanechnikov: kernel_function_ = &epanechnikov_kernel; break;
        default: kernel_function_ = &uniformKernel; break;
        }
    }

    void KDE::kernel_density_estimate() {
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

    void KDE::feed_data(const detail::CircularBuffer &data) {
        std::copy(data.begin(), data.end(), data_.begin());

#ifdef PER4M_USE_FFTW
        if (bandwidth_type_ == BandwidthType::silverman) {
            bandwidth_ = Bandwidth::silverman_1d(data_);
        } else {
            try {
                bandwidth_ = Bandwidth::isj_1d(data_);
            } catch (const std::exception &e) {
                std::cerr << "ISJ bandwidth selection failed: '" << e.what() << "'. Falling back to Silverman's rule." << std::endl;
                bandwidth_ = Bandwidth::silverman_1d(data_);
            }
        }
#else
        bandwidth_ = Bandwidth::silverman_1d(data_);
#endif
    }

    double KDE::estimate(const double x_value) {
        double current_val = lb_;
        const double delta_val = (x_value - lb_) / n_queries_;

        for (int i = 0; i < n_queries_; ++i) {
            x_values_[i] = current_val;
            current_val += delta_val;
        }

        kernel_density_estimate();

        double cdf_value = 0.0;
        for (uint i = 1; i < n_queries_; ++i) {
            cdf_value += 0.5 * (pdf_values_[i] + pdf_values_[i - 1]) * (x_values_[i] - x_values_[i - 1]);
        }
        return cdf_value;
    }

    Kernel get_kernel(const std::string_view kernel) {
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

    BandwidthType get_bandwidth_type(const std::string_view bandwidth_type) {
        if (icase_compare(bandwidth_type, "silverman")) {
            return BandwidthType::silverman;
        }

#ifdef PER4M_USE_FFTW
        if (icase_compare(bandwidth_type, "isj")) {
            return BandwidthType::isj;
        }
#endif

        throw std::invalid_argument("Bandwidth type not defined");
    }
} // namespace per4m