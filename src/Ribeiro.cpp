//
// Created by daniele on 10/02/26.
//

#include <per4m/Ribeiro.h>

#include <cmath>

namespace ribeiro {
    Ribeiro::Ribeiro(double threshold, double improve_pct, double lb, double ub, size_type iterations)
        : threshold_(threshold), improve_pct_(improve_pct), lb_(lb), ub_(ub), iterations_(iterations), data_(iterations) {}

    void Ribeiro::add(double cost) {
#ifdef PER4M_THREAD_SAFE
        std::lock_guard<std::mutex> _(mtx_);
#endif

        data_.add(cost);
        min_ = std::min(min_, cost);
        if (data_.full()) {
            stop_ = get_probability() < threshold_;
        }
    }

    double Ribeiro::phi(const double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }

    std::pair<double, double> Ribeiro::get_mean_and_std() const {
        double sum = 0.0;
        size_type n = 0;

        for (auto x : data_) {
            sum += x;
            ++n;
        }

        double mean = sum / n;
        double variance = 0.0;

        for (const auto x : data_) {
            variance += std::pow(x - mean, 2);
        }
        variance /= n;
        return {mean, std::sqrt(variance)};
    }

    double Ribeiro::get_probability() const {
        if (min_ <= lb_) {
            return 0.0;
        }
        if (min_ >= ub_) {
            return 1.0;
        }

        auto [m, s] = get_mean_and_std();

        const double z_lb = (lb_ - m) / s;
        const double z_ub = (ub_ - m) / s;
        const double z_min = (min_ - m) / s;

        const auto phi_z_ub = ub_ == std::numeric_limits<double>::max() ? 1.0 : phi(z_ub);

        const double denom = phi_z_ub - phi(z_lb);
        if (denom <= 0.0)
            return 0.0; // caso patologico

        return (phi(z_min) - phi(z_lb)) / denom;
    }

    bool Ribeiro::stop() const {
#ifdef PER4M_THREAD_SAFE
        std::lock_guard<std::mutex> _(mtx_);
#endif
        return stop_;
    }
} // namespace ribeiro