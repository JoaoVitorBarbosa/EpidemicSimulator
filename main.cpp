/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: joao
 *
 * Created on 8 de Maio de 2017, 18:29
 */

// git: https://bitbucket.org/JV1992/epidemicsimulator

#include <iostream>
#include "RandomGenerator.h"
#include "Vertex.h"
#include <iomanip>
#include <limits>
#include "EpidemicAnalysis.h"
#include <map>
#include "EpidemicManager.h"
#include "Params.h"
#include "json.hpp"
#include <thread>
#include "Logger.h"

using namespace std;

// json
// https://github.com/nlohmann/json
using json = nlohmann::json;

// static properties must be defined after class declaration
LogLevel Logger::level = LogLevel::Debug;

void to_json(json& j, const Params& p) {
    j = json{
        {"Runs", p.Runs},
        {"Time", p.Time},
        {"OutputDir", p.OutputDir},
        {"Operation", p.Operation},
        {"p", p.p},
        {"KS", p.KS},
        {"TimeS", p.TimeS},
        {"DebugLevel", p.DebugLevel},
        {"GraphParam",
            {
                {"Type", p.Graph.Type},
                {"n", p.Graph.n},
                {"n2", p.Graph.n2},
                {"path", p.Graph.path},
                {"NS", p.Graph.NS}
            }},
        {"RwParam",
            {
                {"k", p.Rw.k},
                {"ki", p.Rw.ki},
                {"kiPer", p.Rw.kiPer},
                {"lambda", p.Rw.lambda},
                {"gama", p.Rw.gama},
                {"tau", p.Rw.tau},
                {"rwParamArray", p.Rw.rwParamArray},
                {"lambdaS", p.Rw.lambdaS}
            }},
        {"VertexParam",
            {
                {"p", p.Vertex.p},
                {"vertexParamArray", p.Vertex.vertexParamArray}
            }}
    };
}

void from_json(const json& j, Params& p) {
    p.Runs = j.at("Runs").get<int>();
    p.Time = j.at("Time").get<int>();
    p.p = j.at("p").get<double>();
    p.OutputDir = j.at("OutputDir").get<std::string>();
    p.Graph.Type = j.at("GraphParam").at("Type").get<int>();
    p.Graph.n = j.at("GraphParam").at("n").get<int>();
    p.Graph.n2 = j.at("GraphParam").at("n2").get<int>();
    p.Graph.path = j.at("GraphParam").at("path").get<std::string>();
    p.Graph.NS = j.at("GraphParam").at("NS").get<std::string>();
    p.Rw.k = j.at("RwParam").at("k").get<int>();
    p.Rw.ki = j.at("RwParam").at("ki").get<int>();
    p.Rw.kiPer = j.at("RwParam").at("kiPer").get<double>();
    p.Rw.lambda = j.at("RwParam").at("lambda").get<double>();
    p.Rw.gama = j.at("RwParam").at("gama").get<double>();
    p.Rw.tau = j.at("RwParam").at("tau").get<double>();
    p.Rw.lambdaS = j.at("RwParam").at("lambdaS").get<std::string>();
    p.Rw.rwParamArray = j.at("RwParam").at("rwParamArray").get<std::string>();
    p.Vertex.p = j.at("VertexParam").at("p").get<double>();
    p.Vertex.vertexParamArray = j.at("VertexParam").at("vertexParamArray").get<std::string>();
    p.KS = j.at("KS").get<std::string>();
    p.TimeS = j.at("TimeS").get<std::string>();
    p.Operation = j.at("Operation").get<std::string>();
    p.DebugLevel = j.at("DebugLevel").get<int>();
}

void validation() {
    double p = 0.37;
    int count = 0;
    int N = 100000;

    RandomGenerator rg;

    for (int i = 0; i < N; i++) {
        if (rg.bernoulli(p))
            count++;
    }

    std::cout << "Count, N: " << count << "," << N << " " << count / (double) N << "%" << std::endl;

}

void test() {
    RandomGenerator rg;
    std::cout << rg.exponential(.8) << std::endl;

    Vertex v1;
    RandomWalk *rw = new RandomWalk();
    rw->set_gama(2.0);

    std::cout << v1.getRandomWalkList().size() << std::endl;

    rw->set_rw_position(v1.setRandomWalk(rw));

    std::cout << v1.getRandomWalkList().size() << std::endl;

    v1.getRandomWalkList().erase(rw->get_rw_position());

    std::cout << v1.getRandomWalkList().size() << std::endl;

    std::cout << "uniform: " << rg.uniform(0, 10) << std::endl;

    delete rw;
}

void testBipartiteGraph() {
    auto graph = GraphGenerator::Bipartite(1, 5);
    std::cout << "Vértices: " << graph.num_vertices << std::endl;
    std::cout << "Arestas: " << graph.num_arestas << std::endl;
    for (int i = 0; i < graph.vetorAdj.size(); i++)
        for (int j = 0; j < graph.vetorAdj[i].size(); j++)
            std::cout << i + 1 << "-" << graph.vetorAdj[i][j] << std::endl;
}

void ploting_expotential(){
    Gnuplot gp;

    gp << "set title 'CCDF tempo gasto no estado \n";
    gp << "set xlabel 'Tempo gasto no estado' \n";
    gp << "set ylabel 'Fração' \n";
    gp << "set term png \n";
    gp << "set output 'teste.png' \n";
    //gp << "set label 'yield point' at 0.0, 0.958 \n";
    //gp << "set yrange [:1] \n";
    gp << "set logscale y \n";
    gp << "plot [0:] exp(-2*x) \n"; // cdf exponential

    //gp.send1d(susceptiblesPts);
}

void validate_torus(){
    int n = 4;
    ManipulaGrafoV graph;
    graph = GraphGenerator::Torus(n);
    for(int i = 1; i <= graph.num_vertices; i++)
    {
        std::cout << i << " | ";
        for(int j = 0; j < graph.get_vizinhos(i).size(); j++)
            std::cout << graph.get_vizinhos(i)[j] << " ";
        std::cout << std::endl;
    }
}

void validate_star(){
    int n1 = 1;
    int n = 10;
    ManipulaGrafoV graph;
    graph = GraphGenerator::Bipartite(n, n1);
    for(int i = 1; i <= graph.num_vertices; i++)
    {
        std::cout << i << " | ";
        for(int j = 0; j < graph.get_vizinhos(i).size(); j++)
            std::cout << graph.get_vizinhos(i)[j] << " ";
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Missing filepath parameter." << std::endl;
        //return 0;
    }

    try {
        //argv[1] = "/home/joao/Mestrado/Simulador/Source\ Code/Epidemic_Simulator/data/Params.json";

        string paramsPath = argv[1];
        std::cout << "Parsing file " << paramsPath << std::endl;
        json j;
        std::ifstream i(paramsPath);
        i >> j;

        // converts json to object
        Params params = j;

        json ar = json::parse(params.Rw.rwParamArray);
        std::vector<std::vector<double>> rws = ar;
        for (int i = 0; i < rws.size(); i++)
            params.Rw.rwParamVector[rws[i][0]] = rws[i];

        json arVertex = json::parse(params.Vertex.vertexParamArray);
        std::vector<std::vector<double>> vertices = arVertex;
        for (int i = 0; i < vertices.size(); i++)
            params.Vertex.vertexParamVector[vertices[i][0]] = vertices[i][1];
        
        json arKS = json::parse(params.KS);
        std::vector<int> KS = arKS;
        params.KS_vector = KS;
        
        json arTimeS = json::parse(params.TimeS);
        std::vector<int> TimeS = arTimeS;
        params.Time_vector = TimeS;
        
        json arNS = json::parse(params.Graph.NS);
        std::vector<int> NS = arNS;
        params.Graph.NS_vector = NS;
        
        json arlambdaS = json::parse(params.Rw.lambdaS);
        std::vector<double> lambdaS = arlambdaS;
        params.Rw.lambda_vector = lambdaS;

        // use parameters of command line if passed
        if(argc > 2)
            params.Time = atof(argv[2]);
        if(argc > 3)
        {
            params.Graph.n = atoi(argv[3]);
            params.Graph.n2 = atoi(argv[3]);
        }
        if(argc > 4)
            params.Rw.k = atoi(argv[4]);
        if(argc > 5)
        {
            params.Rw.kiPer = atof(argv[5]);
        }
        if(argc > 6)
            params.Runs = atoi(argv[6]);
                           
        std::string op = "sim";
        op = params.Operation;
        
        if(argc > 7)
            op = argv[7];
        
        if(argc > 8)
            params.OutputDir = argv[8];
        
        if(argc > 9)
            params.Rw.lambda = atof(argv[9]);
        
        if(argc > 10)
            params.Rw.gama = atof(argv[10]);
        
        
        params.Rw.ki = params.Rw.k * params.Rw.kiPer;
        
        j = params;
        
        Logger::level = (LogLevel) params.DebugLevel;
        
        EpidemicManager manager(false);
        std::string jsonDump = j.dump(4);
        
        if(op == "sim")
            manager.start_simulation(params, jsonDump, true);
        else if(op == "sim_by_k")
            manager.mean_time_epidemic_by_k(params, jsonDump);
        else if(op == "sim_by_n")
            manager.mean_time_epidemic_by_n(params, jsonDump);
        else if(op == "sim_by_lambda")
            manager.mean_time_epidemic_by_lambda(params, jsonDump);
        else if(op == "scalability_k")
            manager.measure_scalability_by_random_walks(params, jsonDump);
        else if(op == "scalability_vertices")
            manager.measure_scalability_by_vertices(params, jsonDump);
        else if(op == "scalability_time")
            manager.measure_scalability_by_time(params, jsonDump);

    } catch (exception& e) {
        std::string str = "Error while parsing parameter file: ";
        str += argv[1];
        str += " Erro: ";
        str += e.what();

        Logger::Error(str);
    }

    return 0;
}

