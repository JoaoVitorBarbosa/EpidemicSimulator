/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "EpidemicManager.h"

EpidemicManager::EpidemicManager() {

}

void EpidemicManager::startSimulation(Params params, std::string strParams) {
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

    time_t now = std::time(0);
    tm *ltm = localtime(&now);
    //std::string fName = std::to_string(ltm->tm_mday) + "-" + std::to_string(1 + ltm->tm_mon) + "-" + std::to_string(1900 + ltm->tm_year) + " " + std::to_string(ltm->tm_hour) + ":" + std::to_string(ltm->tm_min) + ":" + std::to_string(1 + ltm->tm_sec);
    std::string fName = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday) + std::to_string(ltm->tm_hour) + std::to_string(ltm->tm_min) + std::to_string(1 + ltm->tm_sec);

    params.OutputDir = params.OutputDir + "/" + fName;
    this->createDirectory(params.OutputDir);

    std::cout << "Starting simulation" << "\n";
    
    for(int i = 0; i < params.Runs; i++)
    {   
        std::cout << "--------- ROUND " << i << " --------- \n";
        
        Params paramsR = params;
        paramsR.OutputDir = paramsR.OutputDir + "/Round " + std::to_string(i);
        this->createDirectory(paramsR.OutputDir);
        
        Simulator * sim = new Simulator(paramsR, strParams, graph);
        sim->process();
        
        timeAccumulator(sim->time);
        infectionsAccumulator(sim->infectionTimes);
        
        delete sim;
        
        std::cout << "--------- END ROUND " << i << " --------- \n";
    }
    
    std::cout << "Tempo médio de epidemia: " << boost::accumulators::mean(timeAccumulator) << std::endl;
    std::cout << "Desvio padrão do tempo de epidemia: " << std::sqrt(boost::accumulators::variance(timeAccumulator)) << std::endl;
    
    std::cout << "Média de infecções: " << boost::accumulators::mean(infectionsAccumulator) << std::endl;
    std::cout << "Desvio padrão do número de infecções: " << std::sqrt(boost::accumulators::variance(infectionsAccumulator)) << std::endl;
    
    /*std::string analysisDir = sim.outputDir + "/Analysis";

    boost::filesystem::path dir(analysisDir.c_str());
    if (boost::filesystem::create_directory(dir)) {
        std::cout << "Analysis directory created successfully" << "\n";
    }

    EpidemicAnalysis ep(analysisDir);
    ep.params = params;
    ep.infectedsGraphic(sim.fileNameInfectInterval);
    ep.randomWalkStateTimeSeries(sim.fileNameNumberRandomWalkStates);

    ep.analysisStateTime(sim.randomWalks.at(0)->fileNameParmResult);*/

}

void EpidemicManager::createDirectory(std::string directory)
{
    boost::filesystem::path dir(directory.c_str());
    if (boost::filesystem::create_directory(dir)) {
        std::cout << "Directory " << directory << " created successfully" << "\n";
    }
}
