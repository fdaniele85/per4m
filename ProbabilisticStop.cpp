#include "ProbabilisticStop.h"

#include <cctype>
#include <algorithm>

namespace ffp {
    ProbabilisticStop::ProbabilisticStop(double threshold, double improve_pct, Kernel kernel_type, double lb,
                                         int iterations, int number_of_queries)
            : threshold_(threshold),
              improve_pct_(improve_pct),
              estimator_type_(kernel_type),
              iterations_(iterations),
              number_of_queries_(number_of_queries),
              estimator_(estimator_type_, lb) {}

    bool ProbabilisticStop::stop() {
        std::lock_guard<std::mutex> _(mtx_);
        if (!fed_) {
            return false;
        }

        double target_val = min_ - (improve_pct_ * min_);
        return estimate(target_val) < threshold_;
    }

    double ProbabilisticStop::estimate(double num_value) const { return estimator_.cdf(num_value, number_of_queries_); }

    void ProbabilisticStop::add(double cost) {
        std::lock_guard<std::mutex> _(mtx_);

        data_.push_back(cost);
        min_ = std::min(min_, cost);

        if (data_.size() % iterations_ == 0) {
            estimator_.feed_data(data_);
            fed_ = true;
        }
    }
} // namespace ffp