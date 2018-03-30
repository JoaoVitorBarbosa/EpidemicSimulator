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
    double kiPer;
    double lambda;
    std::string lambdaS;
    std::vector<double> lambda_vector;
    double gama;
    double tau;
    std::string rwParamArray; // gets array of json from parameters file
    std::map<int, std::vector<double>> rwParamVector; //for each rw, store its parameters [{rw, [code, vertex, lambda, gama, tau, state (0-Sus, 1-Cont, 2-Inf)]}]
};

struct VertexParam{
    double p;
    std::string vertexParamArray; // gets array of json from parameters file
    std::map<int, double> vertexParamVector; // stores vertex, parameter p.
};

enum GraphType{
    File = 0,
    Ring = 1,
    Torus = 2,
    Clique = 3, 
    Bipartite = 4 // to create a star, pass n = 1 and n2 >= 1
};

struct GraphParam {
    int Type; // 0-File, 1-Ring, 2-Torus, 3-Clique, 4-Bipartite
    int n;
    int n2; // parameter for bipartite graph
    std::string path;
    std::string NS;
    std::vector<int> NS_vector;
};


/// Struct used to convert json parameters in c++ object 
struct Params {
    int Time;
    int Runs;
    double p;
    int DebugLevel;
    bool cont_simulation;
    GraphParam Graph;
    RwParam Rw;
    VertexParam Vertex;
    std::string OutputDir;
    std::string KS;
    std::string TimeS;
    std::vector<int> KS_vector;
    std::vector<int> Time_vector;
    std::string Operation;
};

#endif /* PARAMS_H */
