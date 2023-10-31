//
// Created by daniele on 31/10/23.
//

#include "../include/ProbabilisticStop.h"

#include <limits>
#include <memory>
#include <mlpack/methods/kde.hpp>

namespace dferone {
    ProbabilisticStop::ProbabilisticStop(double threshold, int sample, double bandwith, double LB, double UB, bool repeated, bool minimizing)
        : threshold_(threshold),
          sample_(sample),
          bandwith_(bandwith),
          lb_(LB),
          ub_(UB),
          repeated_(repeated),
          minimizing_(minimizing),
          incumbent_(minimizing ? std::numeric_limits<double>::max() : std::numeric_limits<double>::min()),
          query_(1, number_of_queries_, arma::fill::zeros) {
        auto step = (ub_ - lb_) / number_of_queries_;

        for(uint i = 0; i < number_of_queries_; ++i) {
            query_(0, i) = lb_ + i * step;
        }
    }

    void ProbabilisticStop::add(double cost) {
        data_.push_back(cost);

        if(minimizing_) {
            incumbent_ = std::min(incumbent_, cost);
        } else {
            incumbent_ = std::max(incumbent_, cost);
        }

        if(data_.size() % sample_ == 0) {
            if (!estimated_ || repeated_) {
                estimations_ = estimateDensity();
                estimated_ = true;
            }
        }
    }
    arma::vec ProbabilisticStop::estimateDensity() {
        arma::mat references(1, data_.size(), arma::fill::zeros);

        for(uint i = 0; i < data_.size(); ++i) {
            references(0, i) = data_[i];
        }

        auto kde = std::make_unique<mlpack::KDEModel>(bandwith_);
        arma::vec estimations;

        // TODO Il kernel dovrebbe essere un parametro
        kde->KernelType() = mlpack::KDEModel::GAUSSIAN_KERNEL;
        mlpack::util::Timers timers;
        kde->BuildModel(timers, std::move(references));

        arma::mat copy = query_;
        kde->Evaluate(timers, std::move(copy), estimations);

        return estimations;
    }

    double ProbabilisticStop::estimateProbability() const {
        return estimateProbability(incumbent_);
    }

    double ProbabilisticStop::estimateProbability(double incumbent) const {
        if (!estimated_) {
            return 1;
        }

        int index = 0;
        while (query_(0, index) <= incumbent) {
            ++index;
        }

        arma::mat to_evaluate(1, index, arma::fill::zeros);
        arma::vec estimations_to_evaluate(index, arma::fill::zeros);

        for (int i = 0; i < index; ++i) {
            to_evaluate(0, i) = query_(0, i);
            estimations_to_evaluate(i) = estimations_(i);
        }

        arma::mat Z = trapz(to_evaluate, estimations_to_evaluate);

        if (minimizing_) {
            return Z(0, 0);
        }

        return 1 - Z(0, 0);
    }
    bool ProbabilisticStop::stop() const {
        if (!estimated_) {
            return false;
        }

        return  estimateProbability() <= threshold_;
    }
} // namespace dferone