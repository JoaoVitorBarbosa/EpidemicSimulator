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

class RandomGenerator {
public:
    double exponential(double);
    bool bernoulli(double);
    int uniform(int, int);
    double uniform(double, double);
};

#endif /* RANDOMGENERATOR_H */

