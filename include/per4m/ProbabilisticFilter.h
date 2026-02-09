//
// Created by Tommaso on 02/04/2024.
//

#pragma once

#include "KDE.h"
#include <mutex>

namespace ffp {
    class ProbabilisticFilter {
    public:
        /**
         * \brief Constructs a ProbabilisticFilter object with the specified parameters.
         *
         * This constructor initializes the filter with the given threshold, kernel type,
         * number of iterations, and number of queries.
         *
         * \param threshold The threshold value used in the filter.
         * \param kernel_type The type of kernel estimator used in the filter.
         * \param iterations The number of iterations for the filter process.
         * \param number_of_queries The number of queries to estimate the KDE.
         */
        ProbabilisticFilter(double threshold, Kernel kernel_type, int iterations, int number_of_queries);

        /**
         * \brief Performs least squares estimation based on the target gap.
         *
         * This function performs the least squares estimation using the provided target gap.
         *
         * \param target_gap The target gap value for the least squares estimation.
         * \return True if the estimation is successful, false otherwise.
         */
        [[nodiscard]] bool perform_ls(double target_gap);

        /**
         * \brief Adds a new gap value to the filter.
         *
         * This function adds a new gap value to the filter for further processing.
         *
         * \param gap The gap value to be added.
         */
        void add(double gap);

    private:
        double threshold_;      ///< The threshold value used in the filter.
        Kernel estimator_type_; ///< The type of kernel estimator used in the filter.
        int iterations_;        ///< The number of iterations for the filter process.
        KDE estimator_;         ///< The kernel density estimator used in the filter.
        CircularBuffer data_;   ///< The circular buffer for storing gap values.
        std::mutex mtx_;        ///< Mutex for thread-safe operations.
        bool fed_{false};       ///< Flag indicating whether the filter has been fed with data.
    };
} // namespace ffp