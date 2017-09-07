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

    int num_Infecteds;
    int num_Inf_Events;
    int timeLimit;
    int rounds;
    double timeLastNumberInfect;
    std::map<int, double> infectedInterval;

    void writeInfectInterval();
    void changeNumberInfected(int num_Infec);
    void processWalk(Event evt);
    void processInfect(Event evt);
    void processHeal(Event evt);
    Event getEventTop();
    void fillVertices(Params params);
    void setupRandomWalks(RwParam rwParams);
    void fillEvents();
    void infect(Vertex * v, Event evt);
    void beInfected(Vertex * v, Event evt);
    std::string eventToString(EventType evt);
    void writeNumberRwStatePerTime();

public:
    double time;
    // Number of infected Random Walks
    int k;
    int infectionTimes;
    std::string fileNameInfectInterval;
    std::string fileNameNumberRandomWalkStates;
    std::string outputDir;
    std::vector<RandomWalk*> randomWalks;
    std::priority_queue<Event> events;
    std::vector<Vertex*> vertices;
    ManipulaGrafoV graph;

    Simulator(int _k, int _timeLim, int _round, std::string _graph);
    Simulator(Params params, std::string jsonStr, ManipulaGrafoV graph);
    ~Simulator();

    void initialize(Params params);
    
    /// Starts Simulation
    void process();
};

#endif /* SIMULATOR_H */

