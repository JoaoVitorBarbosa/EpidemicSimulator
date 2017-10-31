/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "EpidemicManager.h"

EpidemicManager::EpidemicManager() {
    aggregate_time_ms = 0;
    epidemic_time = 0.0;
    do_analysis = false;
}

EpidemicManager::EpidemicManager(bool _do_analysis){
    aggregate_time_ms = 0;
    epidemic_time = 0.0;
    do_analysis = _do_analysis;    
}

ManipulaGrafoV EpidemicManager::create_graph(Params params) {
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

    return graph;
}

void EpidemicManager::start_simulation(Params params, std::string strParams) {
    ManipulaGrafoV graph = create_graph(params);

    std::string fName = Utils::datetime_now_to_string();

    params.OutputDir = params.OutputDir + "/" + fName;
    Utils::createDirectory(params.OutputDir);

    std::cout << "Starting simulation" << "\n";

    int num_threads = params.Runs;

    //std::thread * threads = new std::thread[num_threads];

    for (int i = 0; i < params.Runs; i++) {
        //threads[i] = std::thread(&EpidemicManager::run_simulation_in_thread, this, params, strParams, graph, i);
        run_simulation_in_thread(params, strParams, graph, i);
    }

    //for (int i = 0; i < params.Runs; i++)
    //        threads[i].join();

    //std::cout << "Duração média da epidemia " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << std::endl;
    std::cout << "N, " << params.Graph.n << " Duração média da epidemia " << epidemic_time / (double) params.Runs << "s" << std::endl;

    //delete [] threads;

    /*std::cout << "Tempo médio de epidemia: " << boost::accumulators::mean(timeAccumulator) << std::endl;
    std::cout << "Desvio padrão do tempo de epidemia: " << std::sqrt(boost::accumulators::variance(timeAccumulator)) << std::endl;

    std::cout << "Média de infecções: " << boost::accumulators::mean(infectionsAccumulator) << std::endl;
    std::cout << "Desvio padrão do número de infecções: " << std::sqrt(boost::accumulators::variance(infectionsAccumulator)) << std::endl;*/

}

void EpidemicManager::run_simulation_in_thread(Params params, std::string strParams, ManipulaGrafoV graph, int run) {
    std::cout << "--------- ROUND " << run << " --------- \n";

    // creating directory
    Params paramsR = params;
    paramsR.OutputDir = paramsR.OutputDir + "/Round " + std::to_string(run);
    Utils::createDirectory(paramsR.OutputDir);

    // get initial time
    std::chrono::high_resolution_clock::time_point begin_time = std::chrono::high_resolution_clock::now();

    // starting simulation
    Simulator * sim = new Simulator(paramsR, strParams, graph);
    sim->process();

    // get final time
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

    // calculate duration of epidemic
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();

    std::mutex mtx;
    mtx.lock();
    aggregate_time_ms += duration;
    epidemic_time += sim->time;
    mtx.unlock();

    std::cout << "Duração da epidemia " << duration / 1000. << "s" << std::endl;

    // starting analysis
    if (do_analysis) {
        timeAccumulator(sim->time);
        infectionsAccumulator(sim->infectionTimes);

        std::string analysisDir = sim->outputDir + "/Analysis";
        Utils::createDirectory(analysisDir);

        EpidemicAnalysis * ep = new EpidemicAnalysis(analysisDir);
        ep->params = params;
        ep->stateDistribution(sim->fileNameInfectInterval, sim->fileNameContractedInterval, sim->fileNameSusceptibleInterval, "System");
        ep->randomWalkStateTimeSeries(sim->fileNameNumberRandomWalkStates);

        ep->outputDir = analysisDir + "/RandomWalks";
        Utils::createDirectory(ep->outputDir);
        for (int i = 0; i < sim->randomWalks.size(); i++) {
            ep->getRandomWalkCCDFAndStatistics(sim->randomWalks.at(i)->fileNameParmResult, std::to_string(i), sim->randomWalks.at(i)->gama, sim->randomWalks.at(i)->tau);
            ep->RandomWalkWalkingCCDF(sim->randomWalks.at(i)->fileNameWalkingTimes, std::to_string(i), sim->randomWalks.at(i)->lambda);
        }

        ep->outputDir = analysisDir + "/Vertices";
        Utils::createDirectory(ep->outputDir);

        // CONTRACTED AND SUSCEPTIBLE NOT CREATED FOR VERTICES
        //for (int i = 0; i < sim->vertices.size(); i++)
        //ep->stateDistribution(sim->vertices.at(i)->fileNameTimeResult,  std::to_string(i));

        ep->outputDir = sim->outputDir + "/Analysis";
        ep->getSystemCCDFAndStatistics(sim->outputDir + "/RandomWalks", sim->k);
        
        delete ep;
    }

    delete sim;
    
    std::cout << "--------- END ROUND " << run << " --------- \n";

}
