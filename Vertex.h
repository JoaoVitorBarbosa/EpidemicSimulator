/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Vertex.h
 * Author: joao
 *
 * Created on 8 de Maio de 2017, 19:58
 */

#ifndef VERTEX_H
#define VERTEX_H

#include<list>
#include "RandomWalk.h"
#include <vector>>

class Vertex {
    // number of rw infecteds
    int rwInfecteds;
    // probability to rw be infected
    double p;
    // rw list. Store pointer to rw instead of vector index
    std::list<RandomWalk*> randomWalks;
    std::vector<std::pair<double, int> > timeNumberInfect;
    std::string fileNameTimeResult;
    int code;
public:

    Vertex();
    Vertex(int, double, int);
    ~Vertex();

    void setP(double);
    double getP();

    void setRwInfecteds(int);
    void decreaseRwInfecteds(double time);
    void increaseRwInfecteds(double time);
    int getRwInfecteds();
    void writeTimeInfected(std::pair<double, int> pair);

    // Set RandomWalk and return iterator that point to position
    std::list<RandomWalk*>::iterator setRandomWalk(RandomWalk* rw);
    std::list<RandomWalk*>& getRandomWalkList();

    void eraseRandomWalk(std::list<RandomWalk*>::iterator it);
};

#endif /* VERTEX_H */

