/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   RandomWalk.h
 * Author: joao
 *
 * Created on 8 de Maio de 2017, 19:50
 */
#ifndef RANDOMWALK_H
#define RANDOMWALK_H

#include<list>
#include "Event.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "gnuplot-iostream.h"

enum State {
    Susceptible = 0,
    Contracted = 1,
    Infected = 2
};

struct History {

    History() {
    }

    History(int _rw, int _vertex, double _time, std::string _state, std::string _event, std::string _effect) {
        rw = _rw;
        vertex = _vertex;
        time = _time;
        state = _state;
        event = _event;
        effect = _effect;
        rw_inf = -1;
    }

    History(int _rw, int _vertex, double _time, std::string _state, std::string _event, std::string _effect, int _rw_inf) {
        rw = _rw;
        vertex = _vertex;
        time = _time;
        state = _state;
        event = _event;
        effect = _effect;
        rw_inf = _rw_inf;
    }

    int rw;
    // rw that infected rw
    int rw_inf;
    int vertex;
    double time;
    std::string state;
    std::string event;
    std::string effect;

    std::string getRWInfected() {
        return rw_inf == -1 ? "" : std::to_string(rw_inf);
    }
};

class RandomWalk {
public:

    // vetex where rw is located
    int vertex;
    // code of rw
    int code;
    // position of rw in vetex rw vector
    std::list<RandomWalk*>::iterator rwIt;
    double lambda;
    double gama;
    double tau;
    State state;
    bool infectEvent; 
    std::string fileName;
    std::string fileNameParmResult;
    std::string outputDir;

    // store time intervals that rw spend in infected state
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsWalking;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsInfected;
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsContracted;

    // store timestamp and state when state transition occurs
    std::vector<std::pair<double, int>> timeStampStatesChange;

    RandomWalk();
    RandomWalk(int v, double lambda, double gama, double tau, State, int code, std::string _outputDir);

    void setVertex(int);
    int getVertex();

    void setCode(int);
    int getCode();

    void setState(State);
    State getState();

    void setLambda(double);
    double getLambda();

    void setGama(double);
    double getGama();

    void setTau(double);
    double getTau();

    void setRwPosition(std::list<RandomWalk*>::iterator);
    std::list<RandomWalk*>::iterator getRwPosition();

    bool isInfected();

    std::string stateToString();
    void insertTimeInfected(double _t);
    void insertTimeContracted(double _t);
    void insertTimeWalking(double _t);
    
    void WriteResults();

    void writeEvent(History history);
    void writeFile(History history);

    void setTimeStateChange(double _t, int _s);

    void drawStates();

};


#endif /* RANDOMWALK_H */

