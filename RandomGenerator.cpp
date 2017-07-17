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

void RandomGenerator::validateExponential(double rate, int n)
{
    double tax = rate;
    std::multimap<double,double> map;
    
    for(int i = 0; i < n; i++)
    {
        double v = this->exponential(tax);
        map.insert(std::pair<double,double>(v,v));
    }
    
    double estM = 1/rate;
    double estVa = estM;
    double estMed = std::log(2)/rate;
    double mean, var, median = 0.0;
    
    int count = 0;
    for (std::multimap<double,double>::iterator it=map.begin(); it!=map.end(); ++it)
    {
        if(count == n/2)
            median = (*it).first;
        count++;
        mean += (*it).first;
    }
    
    std::cout << "---------- " << "Statistical " << "Sample   " << std::endl;
    std::cout << "Mean: " << estM << "   |   " << mean/n << std::endl;
    std::cout << "Variance: " << estVa << "   |   " << mean/n << std::endl;
    std::cout << "Median: " << estMed << "   |   " << median << std::endl;
}

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