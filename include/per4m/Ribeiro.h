//
// Created by daniele on 10/02/26.
//

#pragma once

#include <per4m/detail/CircularBuffer.h>
#include <limits>
#include <tuple>
#ifdef PER4M_THREAD_SAFE
#include <mutex>
#endif

namespace ribeiro {

    class Ribeiro {
    public:
        using size_type = per4m::detail::CircularBuffer::size_type;
        /// Constructor
        ///
        /// \param threshold and \param improve_pct  regulate the stopping of the probabilistic stop. If
        /// the probability of improving at least improve_pct is lower than threshold, the stoppping criterion is satisfied.
        /// \param lb lower bound
        /// \param ub upper bound
        /// \param iterations the probability is updated every iterations elements
        Ribeiro(double threshold, double improve_pct, double lb, double ub, size_type iterations);
        [[nodiscard]] bool stop() const;
        void add(double cost);



    private:
        double threshold_;
        double improve_pct_;
        double lb_;
        double ub_;
        size_type iterations_;
        double min_{std::numeric_limits<double>::max()};
        per4m::detail::CircularBuffer data_;

        static double phi(double z);
        [[nodiscard]] std::pair<double, double> get_mean_and_std() const;
        [[nodiscard]] double get_probability() const;
        bool stop_ {false};

#ifdef PER4M_THREAD_SAFE
        mutable std::mutex mtx_;
#endif
    };

} // namespace ribeiro