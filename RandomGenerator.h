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

/// Generate random variables according some distribution
class RandomGenerator {
public:
    
    /// Generate exponential variable
    /// \param lambda Exponential parameter
    /// \return 
    double exponential(double lambda);
    
    /// Validate data generate for this exponential. It compares the statistics mean, variance and median with sample
    /// mean, variance and median;
    /// \param rate
    /// \param n samples
    void validateExponential(double rate, int n = 100000);
    
    /// Return true or false according Bernoulli
    /// \param p Bernoulli parameter
    /// \return Return true or false
    bool bernoulli(double p);
    
    /// Generate discrete uniform variable 
    /// \param a
    /// \param b
    /// \return Returns integer between a and b
    int uniform(int a, int b);
    
    /// Generates continuous Uniform variable
    /// \param a
    /// \param b
    /// \return Return double between a and b
    double uniform(double a, double b);
};

#endif /* RANDOMGENERATOR_H */

