/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   EpidemicAnalysis.h
 * Author: joao
 *
 * Created on 13 de Junho de 2017, 21:03
 */

#ifndef EPIDEMICANALYSIS_H
#define EPIDEMICANALYSIS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <set>
#include "gnuplot-iostream.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "Params.h"


class EpidemicAnalysis{
private:
    
    /// Generates statistics of states distributions for each random walk. It also generates the CCDF graphic.
    /// \param timeStampStatesChange Pairs of time interval and State
    /// \param rw Random walk Id
    void createStateStatisticsAndCCDF(std::vector<std::pair<double, int> > timeStampStatesChange, std::string rw);
    
    /// Given a ordered list of elements, calculates the ccdf of the numbers
    /// \param elements
    /// \return Vector of pairs containing the element and the fraction of less elements in the list
    std::vector<std::pair<double,double> > getCCDFPoints(std::multiset<double> elements);
    
    /// Read file and extract randomwalk's walk intervals 
    /// \param filename
    /// \param timeStampWalking Ordered list of walk time
    void readWalkingTime(std::string filename, std::multiset<double> &timeStampWalking);
    
    /// Extract from file, pairs of state and interval
    /// \param filepath
    /// \return Return vector of pairs of (state, interval)
    std::vector<std::pair<int,double>> getStateInterval(std::string filepath);
    
public:
    EpidemicAnalysis(std::string _outputDir);
    Params params;
    
    // Directory to save generated data
    std::string outputDir;
    
    /// Read csv containing timestamp of each state change
    /// \param filename stores the timestamp changed
    /// \timeStampStatesChange ref of vector to stores time and state
    void readTimestampStateChangeCSV(std::string filename, std::vector<std::pair<double, int> > &timeStampStatesChange);
    
    /// Creates graphic with ccdf distribution of states: susceptible, contracted, infected
    /// \param filepath Stores state change and time
    void getRandomWalkCCDFAndStatistics(std::string filepath, std::string title);

    /// Generate statistics and CCDF of all states for all random walks
    /// \param dir Directory to find random walk files
    /// \param num_rw Total number of random walks
    void getSystemCCDFAndStatistics(std::string dir, int num_rw);
    
    /// Creates a graphic with the number of random walk in each state over time
    /// \param file (csv) stores snapshot of (time, rw_infec, rw_cont, rw_inf)
    void randomWalkStateTimeSeries(std::string file);
    
    /// Creates graphic distribution of infected, contracted and susceptible times
    /// \param filepathInfected
    /// \param filepathContracted
    /// \param filepathSusceptible
    /// \param title
    void stateDistribution(std::string filepathInfected, std::string filepathContracted, std::string filepathSusceptible, std::string title);
    
    /// Generate CCDF graphic for walk time of one random walk
    /// \param filepath
    /// \param randomWalk Id of random walk
    void RandomWalkWalkingCCDF(std::string filepath, std::string randomWalk);
    
    
};


#endif /* EPIDEMICANALYSIS_H */

