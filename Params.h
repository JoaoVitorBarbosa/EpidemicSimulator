/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Params.h
 * Author: joao
 *
 * Created on 9 de Julho de 2017, 22:10
 */

#ifndef PARAMS_H
#define PARAMS_H
#include<map>
#include<vector>

struct RwParamArray{
    double lambda;
    double gama;
    double tau;
    int code;
    int vetex;
};

struct RwParam {
    int k;
    int ki;
    double lambda;
    double gama;
    double tau;
    std::string rwParamArray;
    std::map<int, std::vector<double>> rwParamVector; //for each rw, store its parameters [{rw, [code, vertex, lambda, gama, tau, state (0-Sus, 1-Cont, 2-Inf)]}]
};

enum GraphType{
    File = 0,
    Ring = 1,
    Torus = 2,
    Clique = 3, 
    Bipartite = 4
};

struct GraphParam {
    int Type; // 0-File, 1-Ring, 2-Torus, 3-Clique, 4-Bipartite
    int n;
    int n2; // parameter for bipartite graph
    std::string path;
};

struct Params {
    int Time;
    int Rounds;
    double p;
    GraphParam Graph;
    RwParam Rw;
    std::string OutputDir;
};

#endif /* PARAMS_H */
