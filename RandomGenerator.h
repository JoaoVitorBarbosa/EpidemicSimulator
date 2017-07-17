/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   RandomGenerator.h
 * Author: joao
 *
 * Created on 15 de Maio de 2017, 23:03
 */

#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>
#include <chrono>
#include <thread>
#include <map>
#include <iostream>
#include <cmath>

class RandomGenerator {
public:
    double exponential(double);
    /// Validate data generate for this exponential. It compares the statistics mean, variance and median with sample
    /// mean, variance and median;
    /// \param rate
    /// \param n samples
    void validateExponential(double rate, int n = 100000);
    bool bernoulli(double);
    int uniform(int, int);
    double uniform(double, double);
};

#endif /* RANDOMGENERATOR_H */

