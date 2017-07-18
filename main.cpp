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
#include "Simulator.h"
#include <iomanip>
#include <limits>
#include "Params.h"
#include "json.hpp"
#include "GraphGenerator.h"
#include "EpidemicAnalysis.h"
#include <map>

using namespace std;

// json
// https://github.com/nlohmann/json
using json = nlohmann::json;

void to_json(json& j, const Params& p) {
    j = json{
        {"Rounds", p.Rounds},
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
            }}
    };
}

void from_json(const json& j, Params& p) {
    p.Rounds = j.at("Rounds").get<int>();
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

int main(int argc, char** argv) {

    auto graph = GraphGenerator::Bipartite(1,5);
    std::cout << "Vértices: " << graph.num_vertices << std::endl;
    std::cout << "Arestas: " << graph.num_arestas << std::endl;
    for(int i =0; i < graph.vetorAdj.size(); i++)
        for(int j=0; j < graph.vetorAdj[i].size(); j++)
            std::cout << i+1 << "-" << graph.vetorAdj[i][j] << std::endl;
    
    return 0;
    
    if (argc < 2) {
        std::cerr << "Missing Parameter file argument." << std::endl;
        //return 0;
    }

    try {
        //argv[1] = "/home/joao/Mestrado/Simulador/Source\ Code/Epidemic_Simulator/data/Params.json";
        string paramsPath = argv[1];
        std::cout << "Parsing file " << paramsPath << std::endl;
        json j;
        std::ifstream i(paramsPath);
        i >> j;

        Params params = j;
        std::cout << "Graph Type: " << params.Graph.Type << std::endl;
        std::cout << "Graph path: " << params.Graph.path << std::endl;
        std::cout << "RW k: " << params.Rw.k << std::endl;
        std::cout << "RW k: " << params.Rw.ki << std::endl;
        std::cout << "Time: " << params.Time << std::endl;
        std::cout << "Array: " << params.Rw.rwParamArray << std::endl;

        json ar = json::parse(params.Rw.rwParamArray);
        std::vector<std::vector<double>> rws = ar;
        for (int i = 0; i < rws.size(); i++)
            params.Rw.rwParamVector[rws[i][0]] = rws[i];

        ManipulaGrafoV graph;

        switch (params.Graph.Type) {
            case GraphType::File:
                graph.lerArquivo(params.Graph.path, false);
                break;
            case GraphType::Ring:
                graph = GraphGenerator::Ring(params.Graph.n);
                break;
            case GraphType::Clique:
                graph = GraphGenerator::Clique(params.Graph.n);
                break;
            case GraphType::Bipartite:
                graph = GraphGenerator::Bipartite(params.Graph.n, params.Graph.n2);
            default:
                break;
        }

        Simulator sim(params, j.dump(4), graph);
        sim.process();

        std::string analysisDir = sim.outputDir + "/Analysis";

        boost::filesystem::path dir(analysisDir.c_str());
        if (boost::filesystem::create_directory(dir)) {
            std::cout << "Analysis directory created successfully" << "\n";
        }

        EpidemicAnalysis ep(analysisDir);
        ep.params = params;
        ep.infectedsGraphic(sim.fileNameInfectInterval);
        ep.randomWalkStateTimeSeries(sim.fileNameNumberRandomWalkStates);
        
        ep.analysisStateTime(sim.randomWalks.at(0)->fileNameParmResult);

    } catch (exception& e) {
        std::cerr << "Error while parsing parameter file: . " << argv[1] << " Erro: " << e.what() << std::endl;
    }

    return 0;
}

