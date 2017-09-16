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
    std::string fName = std::to_string(1900 + ltm->tm_year) + "-" + std::to_string(1 + ltm->tm_mon) + "-" + std::to_string(ltm->tm_mday) + "-" + std::to_string(ltm->tm_hour) + std::to_string(ltm->tm_min) + std::to_string(1 + ltm->tm_sec);

    params.OutputDir = params.OutputDir + "/" + fName;
    Utils::createDirectory(params.OutputDir);

    std::cout << "Starting simulation" << "\n";

    int num_threads = params.Runs;
    
    std::thread * threads = new std::thread[num_threads];

    for (int i = 0; i < params.Runs; i++) {
        threads[i] = std::thread(&EpidemicManager::runThreadSimulation, this, params, strParams, graph, i);
    }

    for (int i = 0; i < params.Runs; i++)
        threads[i].join();

    delete [] threads;

    std::cout << "Tempo médio de epidemia: " << boost::accumulators::mean(timeAccumulator) << std::endl;
    std::cout << "Desvio padrão do tempo de epidemia: " << std::sqrt(boost::accumulators::variance(timeAccumulator)) << std::endl;

    std::cout << "Média de infecções: " << boost::accumulators::mean(infectionsAccumulator) << std::endl;
    std::cout << "Desvio padrão do número de infecções: " << std::sqrt(boost::accumulators::variance(infectionsAccumulator)) << std::endl;

}

void EpidemicManager::runThreadSimulation(Params params, std::string strParams, ManipulaGrafoV graph, int run) {
    std::cout << "--------- ROUND " << run << " --------- \n";

    // creating directory
    Params paramsR = params;
    paramsR.OutputDir = paramsR.OutputDir + "/Round " + std::to_string(run);
    Utils::createDirectory(paramsR.OutputDir);

    // starting simulation
    Simulator * sim = new Simulator(paramsR, strParams, graph);
    sim->process();

    // starting analysis
    timeAccumulator(sim->time);
    infectionsAccumulator(sim->infectionTimes);

    std::string analysisDir = sim->outputDir + "/Analysis";
    Utils::createDirectory(analysisDir);

    EpidemicAnalysis * ep = new EpidemicAnalysis(analysisDir);
    ep->params = params;
    ep->infectedGraphic(sim->fileNameInfectInterval, sim->fileNameContractedInterval, sim->fileNameSusceptibleInterval,  "System");
    ep->randomWalkStateTimeSeries(sim->fileNameNumberRandomWalkStates);

    ep->outputDir = analysisDir + "/RandomWalks";
    Utils::createDirectory(ep->outputDir);
    for (int i = 0; i < sim->randomWalks.size(); i++)
    {
        ep->analysisStateTime(sim->randomWalks.at(i)->fileNameParmResult, std::to_string(i));
        ep->RandomWalkWalkingCCDF(sim->randomWalks.at(i)->fileNameWalkingTimes, std::to_string(i));
    }

    ep->outputDir = analysisDir + "/Vertices";
    Utils::createDirectory(ep->outputDir);
    
    // CONTRACTED AND SUSCEPTIBLE NOT CREATED TO VERTICES
    //for (int i = 0; i < sim->vertices.size(); i++)
        //ep->infectedGraphic(sim->vertices.at(i)->fileNameTimeResult,  std::to_string(i));

    ep->outputDir = analysisDir;
    ep->analysisAll(sim->outputDir, sim->k);

    delete sim;
    delete ep;

    std::cout << "--------- END ROUND " << run << " --------- \n";
}
