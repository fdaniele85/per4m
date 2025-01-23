//
// Created by Tommaso on 02/04/2024.
//

#pragma once

#include <mutex>
#include <vector>
#include "KDE.h"

namespace ffp {
    class ProbabilisticStop {
    public:
        /// Constructor
        ///
        /// \param threshold and \param improve_pct  regulate the stopping of the probabilistic stop. If
        /// the probability of improving at least improve_pct is lower than threshold, the stoppping criterion is satisfied.
        /// \param kernel_type specifies the type of kernel to be used in the KDE
        /// \param lb lower bound
        /// \param iterations the probability is updated every iterations elements
        /// \param number_of_queries
        ProbabilisticStop(double threshold, double improve_pct, Kernel kernel_type, double lb, int iterations, int number_of_queries);

        [[nodiscard]] bool stop();
        void add(double cost);

    private:


        double threshold_;
        double improve_pct_;
        Kernel estimator_type_;
        int iterations_;
        int number_of_queries_;

        [[nodiscard]] double estimate(double num_value) const;

        KDE estimator_;

        double min_{std::numeric_limits<double>::max()};

        std::vector<double> data_;
        std::mutex mtx_;
        bool fed_{false};
    };
}