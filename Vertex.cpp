/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include <cstddef>

#include "Vertex.h"


Vertex::Vertex(){}

Vertex::Vertex(int _rwInf, double _p, int _code) : rwInfecteds(_rwInf), p(_p), code(_code)
{
//    fileNameTimeResult = "../data/output/" + std::to_string(code) + "_vertexResults.txt";
//    
//    std::ofstream arq; 
//    arq.open(fileNameTimeResult);
//    arq << "#" << '\n';
//    arq.close();
}

Vertex::~Vertex()
{
    //for(std::list<RandomWalk*>::iterator it = randomWalks.begin(); it != randomWalks.end(); it++)
//        delete *it;
}

void Vertex::setP(double _p)
{
    p = _p;
}

double Vertex::getP()
{
    return p;
}

void Vertex::setRwInfecteds(int rw_i)
{
    rwInfecteds = rw_i;
}

int Vertex::getRwInfecteds(){
    return rwInfecteds;
}

void Vertex::decreaseRwInfecteds(double time)
{
    std::pair<double, int> pair = std::make_pair(time, rwInfecteds);
    timeNumberInfect.push_back(pair);
    writeTimeInfected(pair);
    rwInfecteds--;
}

void Vertex::increaseRwInfecteds(double time)
{
    std::pair<double, int> pair = std::make_pair(time, rwInfecteds);
    timeNumberInfect.push_back(pair);
    writeTimeInfected(pair);
    rwInfecteds++;
}

void Vertex::writeTimeInfected(std::pair<double, int> pair)
{
    /*std::ofstream arq; 
    arq.open(fileNameTimeResult, std::ofstream::out | std::ofstream::app);

    arq << std::fixed << pair.first << "," << pair.second << '\n';
    
    arq.close();*/
}


std::list<RandomWalk*>::iterator Vertex::setRandomWalk(RandomWalk* rw)
{
    randomWalks.push_front(rw);
    std::list<RandomWalk*>::iterator it = randomWalks.begin();
    return it;
}

std::list<RandomWalk*>& Vertex::getRandomWalkList(){
    return randomWalks;
}

void Vertex::eraseRandomWalk(std::list<RandomWalk*>::iterator it)
{
    if(!randomWalks.empty())
    {
        randomWalks.erase(it);
    }
}