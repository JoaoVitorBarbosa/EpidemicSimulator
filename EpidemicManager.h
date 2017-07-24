/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   EpidemicManager.h
 * Author: joao
 *
 * Created on 23 de Julho de 2017, 19:54
 */

#ifndef EPIDEMICMANAGER_H
#define EPIDEMICMANAGER_H

#include "Params.h"
#include "Simulator.h"
#include "GraphGenerator.h"
#include "EpidemicAnalysis.h"
#include <unistd.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <ctime>

using namespace std;


class EpidemicManager{
public:
    Params params;
    
    EpidemicManager();
    
    /// Starts Simulation
    /// \param params parameters for simulation
    /// \param paramsStr parameters in string form to print
    void startSimulation(Params params, std::string paramsStr);
    
private:
    void createDirectory(std::string directory);
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> timeAccumulator;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> infectionsAccumulator;
};

#endif /* EPIDEMICMANAGER_H */

