//
// Created by daniele on 31/10/23.
//

#pragma once

#include <vector>
#include <mlpack/methods/kde/kde_model.hpp>
#include <armadillo>

namespace dferone {
    class ProbabilisticStop {
    public:
        /// Constructor
        ///
        /// \param threshold If estimated probability is lower than threshold, the stoppping criterion is satisfied
        /// \param sample Number of solutions use to estimateDensity probability
        /// \param bandwith Bandiwth for the Kernel Density Estimator
        /// \param LB Lower bound
        /// \param UB Upper bound
        /// \param repeated True if estimation is repeated each sample values, False to estimate just once
        /// \param minimizing True if minimizing problem, false in case of maximum problem
        ProbabilisticStop(double threshold, int sample, double bandwith, double LB, double UB, bool repeated, bool minimizing = true);

        /// Checks the stopping condition
        /// \return True if the estimated probability is lower than threshold, false otherwise
        [[nodiscard]] bool stop() const;

        /// Add a new value to the sample
        /// \param cost Cost of the new solution
        void add(double cost);

        [[nodiscard]] double estimateProbability() const;
        [[nodiscard]] double estimateProbability(double incumbent) const;

    private:
        double threshold_;
        int sample_;
        double bandwith_;
        double lb_;
        double ub_;
        bool repeated_;
        bool minimizing_;
        std::vector<double> data_;
        double incumbent_;
        bool estimated_ { false };
        arma::mat query_;
        arma::vec estimations_;

        // FIXME da aggiungere come parametro?
        constexpr static int number_of_queries_ = 5000;

        arma::vec estimateDensity();
    };

} // namespace dferone