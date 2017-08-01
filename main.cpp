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

using namespace std;

// json
// https://github.com/nlohmann/json
using json = nlohmann::json;

void to_json(json& j, const Params& p) {
    j = json{
        {"Runs", p.Runs},
        {"Time", p.Time},
        {"OutputDir", p.OutputDir},
        {"p", p.p},
        {"GraphParam",
            {
                {"Type", p.Graph.Type},
                {"n", p.Graph.n},
                {"path", p.Graph.path}
            }},
        {"RwParam",
            {
                {"k", p.Rw.k},
                {"ki", p.Rw.ki},
                {"lambda", p.Rw.lambda},
                {"gama", p.Rw.gama},
                {"tau", p.Rw.tau},
                {"rwParamArray", p.Rw.rwParamArray}
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
    p.Graph.path = j.at("GraphParam").at("path").get<std::string>();
    p.Rw.k = j.at("RwParam").at("k").get<int>();
    p.Rw.ki = j.at("RwParam").at("ki").get<int>();
    p.Rw.lambda = j.at("RwParam").at("lambda").get<double>();
    p.Rw.gama = j.at("RwParam").at("gama").get<double>();
    p.Rw.tau = j.at("RwParam").at("tau").get<double>();
    p.Rw.rwParamArray = j.at("RwParam").at("rwParamArray").get<std::string>();
    p.Vertex.p = j.at("VertexParam").at("p").get<double>();
    p.Vertex.vertexParamArray = j.at("VertexParam").at("vertexParamArray").get<std::string>();
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
    rw->setGama(2.0);

    std::cout << v1.getRandomWalkList().size() << std::endl;

    rw->setRwPosition(v1.setRandomWalk(rw));

    std::cout << v1.getRandomWalkList().size() << std::endl;

    v1.getRandomWalkList().erase(rw->getRwPosition());

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


int main(int argc, char** argv) {

    /*std::thread first([]() {
        int i = 0;
        for (int j = 0; j < 1000000; j++) {
            i = j;
        }
    });
    std::thread second([]() {
        int i = 0;
        for (int j = 0; j < 1000000; j++) {
            i = j;
        }
    });
    
    std::thread third([]() {
        int i = 0;
        for (int j = 0; j < 1000000; j++) {
            i = j;
        }
    });

    std::cout << "main, foo and bar now execute concurrently...\n" << std::endl;


    first.join();
    second.join();
    third.join();

    std::cout << "foo and bar completed.\n";

    return 0;*/

    if (argc < 2) {
        std::cerr << "Missing Parameter file argument." << std::endl;
        //return 0;
    }

    try {
        argv[1] = "/home/joao/Mestrado/Simulador/Source\ Code/Epidemic_Simulator/data/Params.json";

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

        EpidemicManager manager;
        manager.startSimulation(params, j.dump(4));

    } catch (exception& e) {
        std::cerr << "Error while parsing parameter file: . " << argv[1] << " Erro: " << e.what() << std::endl;
    }

    return 0;
}

