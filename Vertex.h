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
#include <map>

class Vertex {
    // number of rw infecteds
    int rwInfecteds;
    // probability to rw be infected
    double p;
    // rw list. Store pointer to rw instead of vector index
    std::list<RandomWalk*> randomWalks;
    std::map<int, double> timeNumberInfect;
    int code;
    double timeLastNumberinfect;
    int kMax;
    std::string dirToSave;
    double totalTimeWithInfc;
    
public:
    std::string fileNameTimeResult;
    
    // directory to save all vertex data
    std::string outputDir;
    
    Vertex();
    Vertex(int rwinf_0, double p, int code, int _kMax, std::string _dirToSave);
    ~Vertex();

    void setP(double);
    double getP();

    void setRwInfecteds(int);
    void decreaseRwInfecteds(double time);
    void increaseRwInfecteds(double time);
    int getRwInfecteds();
    
    /// Write file with infected distribution
    /// \param time is the last time of system
    void writeTimeInfected(double time);

    /// Set RandomWalk and return iterator that point to its position
    /// \param rw
    /// \return 
    std::list<RandomWalk*>::iterator setRandomWalk(RandomWalk* rw);
    std::list<RandomWalk*>& getRandomWalkList();

    /// Remove Random Walk from RW list using its iterator
    /// \param it RW Iterator
    void eraseRandomWalk(std::list<RandomWalk*>::iterator it);
};

#endif /* VERTEX_H */

