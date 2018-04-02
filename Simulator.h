/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Simulator.h
 * Author: joao
 *
 * Created on 8 de Maio de 2017, 21:56
 */


#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <math.h>
#include <map>
#include <vector>
#include <queue>
#include "Utils.h"
#include "RandomWalk.h"
#include "ManipulaGrafo.h"
#include "Vertex.h"
#include "RandomGenerator.h"
#include "Event.h"
#include "Params.h"
#include <unistd.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <ctime>


class Simulator {

    RandomGenerator rg;
    int num_Infected;
    int num_Contracted;
    int num_Inf_Events;
    int rounds;
    bool continue_simulation_after_epidemic;
    double time_of_last_number_of_infected;
    /// Stores Time interval that system had k infected
    std::map<int, double> infected_intervals;
    double time_of_last_number_of_contracted;
    /// Stores Time interval that system had k contracted
    std::map<int, double> contracted_intervals;
    double time_of_last_number_of_susceptible;
    /// Stores Time interval that system had k suscetible
    std::map<int, double> susceptible_intervals;

    void write_infected_intervals();
    void write_contracted_intervals();
    void write_susceptible_intervals();
    
    /// Changes number of infected e compute time interval. Increase or decrease by 1.
    /// \param num_Infec New number of infected
    void change_infected_number(int num_Infec);
    
    /// Changes number of contracted e compute time interval simulator spends with this number. Increase or decrease by 1.
    /// \param num_Contracted New number of contracted
    void changeNumberContracted(int num_Contracted);
    
    void process_walk(Event evt);
    void process_infect(Event evt);
    void process_heal(Event evt);
    Event get_top_event();
    void fill_vertices(Params params);
    void setupRandomWalks(RwParam rwParams);
    void fillEvents();
    void infect(Vertex * v, Event evt);
    void beInfected(Vertex * v, Event evt);
    std::string eventToString(EventType evt);
    void writeNumberRwStatePerTime(std::string evt);

    std::list<std::pair<double, double> > infected_density_over_time;
    // store pairs of (time, # infected). It's used to create graph of infected in function of time (and for average over runs)
    std::string file_infected_density;
    void write_infected_density();
    
public:
    int limit_time_epidemic;
    double time;
    std::list<std::pair<double, int> > infected_time;
    // Number of infected Random Walks
    int k;
    int infectionTimes;
    std::string fileNameInfectInterval;
    std::string fileNameContractedInterval;
    std::string fileNameSusceptibleInterval;
    std::string fileNameNumberRandomWalkStates;
    std::string file_name_system_results;
    std::string outputDir;
    std::vector<RandomWalk*> randomWalks;
    std::priority_queue<Event> events;
    std::vector<Vertex*> vertices;
    ManipulaGrafoV graph;
    bool do_analysis;
    
    Simulator(int _k, int _timeLim, int _round, std::string _graph);
    Simulator(Params params, std::string jsonStr, ManipulaGrafoV graph);
    ~Simulator();

    void initialize(Params params);
    
    /// Starts Simulation
    void process();
};

#endif /* SIMULATOR_H */

