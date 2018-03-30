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

#include <list>
#include "RandomWalk.h"
#include <map>

class Vertex {
    // number of rw infecteds
    int rwInfecteds;
    // probability to rw be infected. Bernoulli parameter
    double p;
    // rw list. Store pointer to rw instead of vector index
    std::list<RandomWalk*> randomWalks;
    std::map<int, double> timeNumberInfect;
    int code;
    double timeLastNumberinfect;
    int kMax;
    std::string dirToSave;
    double totalTimeWithInfc;
    
    // validation
    int total_encounters;
    int total_encounters_with_transmission;
    
    // store sum o k (infected rw) to encounters that generate infection 
    int sum_success_k;
    // store sum o k (infected rw) to encounters that don't generate infection 
    int sum_fail_k;
    
    double get_empiric_p_type_1();
    
public:
    std::string fileNameTimeResult;
    
    // directory to save all vertex data
    std::string outputDir;
    bool do_analysis;
    
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
    
    void increase_encounters_by(int n);
    void increase_encounters_with_transmition_by(int n);
    
    void sum_fail_encounters(int k);
    void sum_success_encounters(int k);
};

#endif /* VERTEX_H */

