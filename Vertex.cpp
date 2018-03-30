/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include <cstddef>

#include "Vertex.h"

Vertex::Vertex() {
}

Vertex::Vertex(int _rwInf, double _p, int _code, int _kMax, std::string _dirToSave) : rwInfecteds(_rwInf), p(_p), code(_code), kMax(_kMax), dirToSave(_dirToSave) {
    do_analysis = true;
    fileNameTimeResult = dirToSave + "/" + std::to_string(code) + "_vertex_results.txt";
    totalTimeWithInfc = 0.0;
    timeLastNumberinfect = 0.0;
    total_encounters = 0;
    total_encounters_with_transmission = 0;
    sum_fail_k = 0;
    sum_success_k = 0;
}

Vertex::~Vertex() {
    //for(std::list<RandomWalk*>::iterator it = randomWalks.begin(); it != randomWalks.end(); it++)
    //        delete *it;
}

void Vertex::setP(double _p) {
    p = _p;
}

double Vertex::getP() {
    return p;
}

void Vertex::setRwInfecteds(int rw_i) {
    rwInfecteds = rw_i;
}

int Vertex::getRwInfecteds() {
    return rwInfecteds;
}

void Vertex::decreaseRwInfecteds(double time) {
    double interval = time - this->timeLastNumberinfect;
    double actInterval = timeNumberInfect[rwInfecteds];

    timeNumberInfect[rwInfecteds] = interval + actInterval;
    totalTimeWithInfc += interval;

    timeLastNumberinfect = time;
    rwInfecteds--;
}

void Vertex::increaseRwInfecteds(double time) {
    double interval = time - this->timeLastNumberinfect;
    double actInterval = timeNumberInfect[rwInfecteds];

    // don't count time with k = 0
    // if size == 1, so it was 0 before
    if (randomWalks.size() > 1) {
        timeNumberInfect[rwInfecteds] = interval + actInterval;
        totalTimeWithInfc += interval;
    }

    timeLastNumberinfect = time;
    rwInfecteds++;
}

void Vertex::writeTimeInfected(double time) {
    if (do_analysis) {
        std::ofstream arq;
        arq.open(fileNameTimeResult, std::ofstream::out | std::ofstream::app);

        for (std::map<int, double>::iterator it = timeNumberInfect.begin(); it != timeNumberInfect.end(); it++)
            arq << it->first << "," << it->second / totalTimeWithInfc << std::endl;

        arq << "Total encounters type 2: " << total_encounters << std::endl;
        arq << "Total encounters type 2 with transmission: " << total_encounters_with_transmission << std::endl;
        arq << "Empiric p type 2: " << total_encounters_with_transmission / (double) total_encounters << std::endl;
        arq << "Empiric p type 1: " << get_empiric_p_type_1() << std::endl;
        arq << "\n";

        for (int i = 0; i <= kMax; i++)
            arq << i << "," << timeNumberInfect[i] / totalTimeWithInfc << std::endl;

        arq.close();
    }
}

std::list<RandomWalk*>::iterator Vertex::setRandomWalk(RandomWalk* rw) {
    randomWalks.push_front(rw);
    std::list<RandomWalk*>::iterator it = randomWalks.begin();
    return it;
}

std::list<RandomWalk*>& Vertex::getRandomWalkList() {
    return randomWalks;
}

void Vertex::eraseRandomWalk(std::list<RandomWalk*>::iterator it) {
    if (!randomWalks.empty()) {
        randomWalks.erase(it);
    }
}

void Vertex::increase_encounters_by(int n) {
    this->total_encounters += n;
}

void Vertex::increase_encounters_with_transmition_by(int n) {
    this->total_encounters_with_transmission += n;
}

void Vertex::sum_fail_encounters(int k) {
    sum_fail_k += k;
}

void Vertex::sum_success_encounters(int k) {
    sum_success_k += k;
}

double Vertex::get_empiric_p_type_1() {
    return sum_fail_k / (double) (sum_success_k - sum_fail_k);
}