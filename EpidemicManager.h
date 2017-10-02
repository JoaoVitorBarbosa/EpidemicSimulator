/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   EpidemicManager.h
 * Author: joao
 *
 * Created on 23 de Julho de 2017, 19:54
 * 
 * Manage epidemic, setting parameters and starting simulation
 * 
 */

#ifndef EPIDEMICMANAGER_H
#define EPIDEMICMANAGER_H

#include "Utils.h"
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
#include <thread>
#include <chrono>
#include <mutex>

using namespace std;


class EpidemicManager{
public:
    Params params;
    
    EpidemicManager();
    
    /// Starts Simulation
    /// \param params parameters for simulation
    /// \param paramsStr parameters in string form to print
    void start_simulation(Params params, std::string paramsStr);
    
private:
    
    /// Create graph according to specified in params
    /// \param params
    /// \return 
    ManipulaGrafoV create_graph(Params params);
    
    /// Run a simulation inside of a thread
    /// \param params
    /// \param paramsStr Parameter that will be printed in log
    /// \param graph 
    /// \param run Represents id of run
    void run_simulation_in_thread(Params params, std::string paramsStr, ManipulaGrafoV graph, int run);
    
    /// duration of all runs
    int aggregate_time_ms;
    
    /// Boost accumulator used to store duration of epidemic
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> timeAccumulator;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> infectionsAccumulator;
};

#endif /* EPIDEMICMANAGER_H */

