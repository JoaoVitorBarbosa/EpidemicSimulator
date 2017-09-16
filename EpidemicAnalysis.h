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
    /// 
    /// \param timeStampStatesChange
    /// \param rw
    void createStateStatistics(std::vector<std::pair<double, int> > timeStampStatesChange, std::string rw);
    
    std::vector<std::pair<double,double> > getCCDFPoints(std::multiset<double> elements);

public:
    EpidemicAnalysis(std::string _outputDir);
    Params params;
    std::string outputDir;
    
    /// Read csv containing timestamp of each state change
    /// \param filename stores the timestamp changed
    /// \timeStampStatesChange ref of vector to stores time and state
    void readTimestampStateChangeCSV(std::string filename, std::vector<std::pair<double, int> > &timeStampStatesChange);
    
    /// Creates graphic with ccdf distribution of states: susceptible, contracted, infected
    /// \param filepath Stores state change and time
    void analysisStateTime(std::string filepath, std::string title);

    void printCCDF(std::multiset<double> elements, std::string title);
    void analysisAll(std::string dir, int num_rw);
    
    /// Creates a graphic with the number of random walk in each state over time
    /// \param arquivo (csv) stores snapshot of time, rw_infec, rw_cont, rw_inf
    void randomWalkStateTimeSeries(std::string arquivo);
    
    void infectedGraphic(std::string filepathInfected, std::string filepathContracted, std::string filepathSusceptible, std::string title);
    
    void RandomWalkWalkingCCDF(std::string filepath, std::string randomWalk);
    
    void readWalkingTime(std::string filename, std::multiset<double> &timeStampWalking);
};


#endif /* EPIDEMICANALYSIS_H */

