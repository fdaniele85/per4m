//
// Created by Tommaso on 02/04/2024.
//

#pragma once

#include "KDE.h"
#include <mutex>

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

        /**
         * \brief Checks if the stopping criterion is satisfied.
         *
         * \return True if the stopping criterion is satisfied, false otherwise.
         */
        [[nodiscard]] bool stop();

        /**
         * \brief Adds a new cost value to the data.
         *
         * \param cost The cost value to be added.
         */
        void add(double cost);

    private:
        double threshold_;      ///< The threshold value for the stopping criterion.
        double improve_pct_;    ///< The percentage improvement required for the stopping criterion.
        Kernel estimator_type_; ///< The type of kernel estimator used in the KDE.
        int iterations_;        ///< The number of iterations for updating the probability.

        KDE estimator_; ///< The kernel density estimator used in the stopping criterion.

        double min_{std::numeric_limits<double>::max()}; ///< The minimum cost value observed.

        CircularBuffer data_; ///< Circular buffer to store the cost values.
        std::mutex mtx_;      ///< Mutex for thread-safe operations.
        bool fed_{false};     ///< Flag indicating whether the data has been fed.
    };
} // namespace ffp