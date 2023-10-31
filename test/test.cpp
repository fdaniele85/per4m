#include <gtest/gtest.h>

#include <ProbabilisticStop.h>
#include <random>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    dferone::ProbabilisticStop ps(0.05, 50, 0.01, -20, 20, true, true);

    std::random_device dv;
    std::normal_distribution<double> dis;

    for (uint i = 0; i < 50000; ++i) {
        ps.add(dis(dv));
        ps.stop()
    }

    while ()
}