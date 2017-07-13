/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "RandomGenerator.h"

double RandomGenerator::exponential(double lambda) {
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 mt(seed);

    std::exponential_distribution<double> exp(lambda);

    return exp(mt);
};

bool RandomGenerator::bernoulli(double p) {
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 mt(seed);

    std::uniform_real_distribution<double> unif(0.0, 1.0);
    double t = unif(mt);
    return t < p;
}

int RandomGenerator::uniform(int a, int b) {
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 mt(seed);

    std::uniform_real_distribution<double> unif(a, b);
    return unif(mt);
}

double RandomGenerator::uniform(double a, double b) {
    int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 mt(seed);

    std::uniform_real_distribution<double> unif(a, b);
    return unif(mt);
}