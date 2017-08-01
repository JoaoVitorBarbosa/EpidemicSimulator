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
    std::string fName = std::to_string(ltm->tm_mday) + "-" + std::to_string(1 + ltm->tm_mon) + "-" + std::to_string(1900 + ltm->tm_year) + " " + std::to_string(ltm->tm_hour) + ":" + std::to_string(ltm->tm_min) + ":" + std::to_string(1 + ltm->tm_sec);
    //std::string fName = std::to_string(1900 + ltm->tm_year) + std::to_string(1 + ltm->tm_mon) + std::to_string(ltm->tm_mday) + std::to_string(ltm->tm_hour) + std::to_string(ltm->tm_min) + std::to_string(1 + ltm->tm_sec);

    params.OutputDir = params.OutputDir + "/" + fName;
    this->createDirectory(params.OutputDir);

    std::cout << "Starting simulation" << "\n";

    std::thread * threads = new std::thread[params.Runs];

    for (int i = 0; i < params.Runs; i++) {
        
        threads[i] = std::thread(&EpidemicManager::runThreadSimulation, this, params, strParams, graph, i);
        
        /*threads[i] = std::thread([](EpidemicManager &_this, Params params, std::string strParams, ManipulaGrafoV graph, int i) {
            std::cout << "--------- ROUND " << i << " --------- \n";

            Params paramsR = params;
            paramsR.OutputDir = paramsR.OutputDir + "/Round " + std::to_string(i);
            _this.createDirectory(paramsR.OutputDir);

            Simulator * sim = new Simulator(paramsR, strParams, graph);
            sim->process();

            _this.timeAccumulator(sim->time);
            _this.infectionsAccumulator(sim->infectionTimes);

            std::string analysisDir = sim->outputDir + "/Analysis";
            _this.createDirectory(analysisDir);

            EpidemicAnalysis * ep = new EpidemicAnalysis(analysisDir);
            ep->params = params;
            ep->infectedsGraphic(sim->fileNameInfectInterval, "System");
            ep->randomWalkStateTimeSeries(sim->fileNameNumberRandomWalkStates);

            ep->outputDir = analysisDir + "/RandomWalks";
            _this.createDirectory(ep->outputDir);
            for (int i = 0; i < sim->randomWalks.size(); i++)
                    ep->analysisStateTime(sim->randomWalks.at(i)->fileNameParmResult, std::to_string(i));

                    ep->outputDir = analysisDir + "/Vertices";
                    _this.createDirectory(ep->outputDir);
                for (int i = 0; i < sim->vertices.size(); i++)
                        ep->infectedsGraphic(sim->vertices.at(i)->fileNameTimeResult, std::to_string(i));

                        //ep->outputDir = analysisDir;
                        //ep->analysisAll(analysisDir, sim->k);

                        delete sim;
                        delete ep;

                        std::cout << "--------- END ROUND " << i << " --------- \n";
                }, this, params, strParams, graph, i);*/
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

    Params paramsR = params;
    paramsR.OutputDir = paramsR.OutputDir + "/Round " + std::to_string(run);
    this->createDirectory(paramsR.OutputDir);

    Simulator * sim = new Simulator(paramsR, strParams, graph);
    sim->process();

    timeAccumulator(sim->time);
    infectionsAccumulator(sim->infectionTimes);

    std::string analysisDir = sim->outputDir + "/Analysis";
    this->createDirectory(analysisDir);

    EpidemicAnalysis * ep = new EpidemicAnalysis(analysisDir);
    ep->params = params;
    ep->infectedsGraphic(sim->fileNameInfectInterval, "System");
    ep->randomWalkStateTimeSeries(sim->fileNameNumberRandomWalkStates);

    ep->outputDir = analysisDir + "/RandomWalks";
    this->createDirectory(ep->outputDir);
    for (int i = 0; i < sim->randomWalks.size(); i++)
        ep->analysisStateTime(sim->randomWalks.at(i)->fileNameParmResult, std::to_string(i));

    ep->outputDir = analysisDir + "/Vertices";
    this->createDirectory(ep->outputDir);
    for (int i = 0; i < sim->vertices.size(); i++)
        ep->infectedsGraphic(sim->vertices.at(i)->fileNameTimeResult, std::to_string(i));

    //ep->outputDir = analysisDir;
    //ep->analysisAll(analysisDir, sim->k);

    delete sim;
    delete ep;

    std::cout << "--------- END ROUND " << run << " --------- \n";
}

void EpidemicManager::createDirectory(std::string directory) {
    boost::filesystem::path dir(directory.c_str());
    if (boost::filesystem::create_directory(dir)) {
        std::cout << "Directory " << directory << " created successfully" << "\n";
    }
}
