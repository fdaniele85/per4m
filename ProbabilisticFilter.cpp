#include "ProbabilisticFilter.h"

namespace ffp {
    ProbabilisticFilter::ProbabilisticFilter(const double threshold, const Kernel kernel_type, const int iterations, const int number_of_queries)
        : threshold_(threshold),
          estimator_type_(kernel_type),
          iterations_(iterations),
          estimator_(estimator_type_, iterations, number_of_queries, 0),
          data_(iterations) {}

    void ProbabilisticFilter::add(const double gap) {
        std::lock_guard<std::mutex> _(mtx_);

        data_.add(gap);

        if (data_.full()) {
            estimator_.feed_data(data_);
            fed_ = true;
        }
    }

    bool ProbabilisticFilter::perform_ls(const double target_gap) {
        std::lock_guard<std::mutex> _(mtx_);
        if (!fed_) {
            return true;
        }

        return (1 - estimator_.estimate(target_gap)) >= threshold_;
    }
} // namespace ffp