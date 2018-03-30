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
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <ctime>
#include <thread>
#include <chrono>
#include <mutex>
#include "json.hpp"
#include "Logger.h"

using namespace std;

struct TimeStatistics{
    double std;
    double mean;
    double min;
    double max;
};


class EpidemicManager{
public:
    Params params;
    
    EpidemicManager();
    EpidemicManager(bool do_analysis);
    
    /// Starts Simulation
    /// \param params parameters for simulation
    /// \param paramsStr parameters in string form to print
    void start_simulation(Params params, std::string paramsStr, bool do_parallel);
    
    void mean_time_epidemic_by_k(Params params, std::string strParams);
    void mean_time_epidemic_by_n(Params params, std::string strParams);
    void mean_time_epidemic_by_lambda(Params params, std::string strParams);
    
    void measure_scalability_by_random_walks(Params params, std::string strParams);
    void measure_scalability_by_vertices(Params params, std::string strParams);
    void measure_scalability_by_time(Params params, std::string strParams);
    
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
    
    TimeStatistics time_epidemic_statistics(double mean);
    
    /// duration of all runs
    int aggregate_time_ms;
    double epidemic_time;
    double avg_epidemic_time;
    double std_epidemic_time;
    std::vector<std::tuple<int,double>> means;
    std::vector<double> means_dur_execution;
    int do_analysis;
    std::string output_dir;
    TimeStatistics duration_standard_desviation(std::vector<double> times, double mean);
    
    
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> time_epidemic_acc;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> time_epidemic_execution_acc;
    
    /// Boost accumulator used to store duration of epidemic
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> timeAccumulator;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance>> infectionsAccumulator;
};

#endif /* EPIDEMICMANAGER_H */

