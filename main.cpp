#define MLPACK_PRINT_INFO
#define MLPACK_PRINT_WARN
#include <ProbabilisticStop.h>
#include <mlpack/methods/kde.hpp>
#include <mlpack/methods/kde/kde_model.hpp>
#include <random>

using namespace arma;
using namespace mlpack;
using namespace mlpack::util;
using namespace std;

int main() {
    //    KDEModel* kde= new KDEModel(0.05);
    //    arma::vec estimations;
    //    arma::mat reference(1, 500, arma::fill::randn);
    //
    //    std::default_random_engine generator;
    //    std::poisson_distribution<int> distribution(4.1);
    //
    //    for (int i = 0; i < 500; ++i) {
    //        reference(0, i) = distribution(generator);
    //    }
    //    arma::mat query(1, 3000, arma::fill::zeros);
    //
    //    auto minimum = reference.min();
    //    auto maximum = reference.max();
    //
    //    for (uint i=0;i<3000;++i){
    //        query(0,i) = +0.005*i;
    //    }
    ////    for (uint i = 0; i < 5000; ++i) {
    ////        reference(0, i) *= 1;
    ////    }
    //    std:cout << query << '\n';
    //    kde->KernelType() = KDEModel::GAUSSIAN_KERNEL;
    //    util::Timers timers;
    //    //kde->Bandwidth();
    //    kde->BuildModel(timers, std::move(reference));
    //
    //    arma::mat copy = query;
    //    kde->Evaluate(timers, std::move(query), estimations);
    //
    //    int index = 0;
    //    std::cout << "min: " << minimum << "\n";
    //    std::cout << "max: " << maximum << "\n";
    //    while (copy(0, index) <= 3) {
    //        ++index;
    //    }
    //
    //    arma::mat to_evaluate(1, index, arma::fill::zeros);
    //    arma::vec estimations_to_evaluate(index, arma::fill::zeros);
    //
    //    for (int i = 0; i < index; ++i) {
    //        to_evaluate(0, i) = copy(0, i);
    //        estimations_to_evaluate(i) = estimations(i);
    //    }
    //
    //    //std::cout << estimations_to_evaluate << "\n";
    //
    //
    //    arma::mat Z = trapz(to_evaluate, estimations_to_evaluate);
    //
    //    std::cout << "Z: " << Z(0, 0) << '\n';
    //
    //    return 0;

    dferone::ProbabilisticStop ps(0.0001, 250, 0.05, -20, 20, true, true);

    std::random_device dv;
    std::normal_distribution<double> dis;

    auto min_value = std::numeric_limits<double>::max();

    for(uint i = 0; i < 50000; ++i) {
        auto value = dis(dv);
        min_value = std::min(min_value, value);

        ps.add(value);
        std::cout << i << ": " << min_value << " -- " << ps.estimateProbability(min_value) << " -- " << ps.stop() << '\n';
    }

    return 0;
}
