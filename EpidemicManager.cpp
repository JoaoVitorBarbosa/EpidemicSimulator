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

EpidemicManager::EpidemicManager(bool _do_analysis) {
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
        case GraphType::Torus:
            graph = GraphGenerator::Torus(params.Graph.n);
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

void EpidemicManager::mean_time_epidemic_by_k(Params params, std::string strParams) {
    std::vector<std::pair<double, double> > k_mean_time;
    std::vector<std::pair<double, double> > k_std_time;
    int N = params.Graph.n;

    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_simluation_by_k";
    Utils::createDirectory(params.OutputDir);
    std::string filename_mean_std_time = params.OutputDir + "/time_std_values.txt";

    std::ofstream arq;

    for (int i = 0; i < params.KS_vector.size(); i++) {
        params.Rw.k = params.KS_vector[i];

        // phase transition
        if (avg_epidemic_time >= params.Time) {
            std::cout << "Phase transition achieved, stop varying k. Setting Time limit for avg_time" << std::endl;

            avg_epidemic_time = params.Time;
            std_epidemic_time = 0.0012;
        } else {
            start_simulation(params, strParams);
        }

        k_mean_time.push_back(std::make_pair(params.KS_vector[i], avg_epidemic_time));
        k_std_time.push_back(std::make_pair(params.KS_vector[i], std_epidemic_time));

        arq.open(filename_mean_std_time, std::ofstream::app);
        arq << params.Rw.k << ";" << Utils::replace_dot_to_comma(avg_epidemic_time) << ";" << Utils::replace_dot_to_comma(std_epidemic_time) << "\n";
        arq.close();
    }

    EpidemicAnalysis * ep = new EpidemicAnalysis(params.OutputDir);
    ep->params = params;
    ep->mean_time_of_epidemic_over_k(k_mean_time, N);
    ep->std_time_of_epidemic_over_k(k_std_time, N);
    delete ep;

}

void EpidemicManager::mean_time_epidemic_by_n(Params params, std::string strParams) {
    std::vector<std::pair<double, double> > n_mean_time;
    std::vector<std::pair<double, double> > n_std_time;
    int K = params.Rw.k;
    std::string outputdir;
    for (int i = 0; i < params.Graph.NS_vector.size(); i++) {
        params.Graph.n = params.Graph.NS_vector[i];
        start_simulation(params, strParams);
        n_mean_time.push_back(std::make_pair(params.Graph.NS_vector[i], avg_epidemic_time));
        n_std_time.push_back(std::make_pair(params.Graph.NS_vector[i], std_epidemic_time));

        // phase transition
        // epidemia não morre para n pequenos
        //        if(avg_epidemic_time >= params.Time)
        //        {
        //            std::cout << "Phase transition achieved, stop varying n" << std::endl;
        //            break;
        //        }
    }

    outputdir = output_dir;
    EpidemicAnalysis * ep = new EpidemicAnalysis(outputdir);
    ep->params = params;
    ep->mean_time_of_epidemic_over_n(n_mean_time, K);
    ep->std_time_of_epidemic_over_n(n_std_time, K);
    delete ep;

}

void EpidemicManager::mean_time_epidemic_by_lambda(Params params, std::string strParams) {
    std::vector<std::pair<double, double> > lambda_mean_time;
    std::vector<std::pair<double, double> > lambda_std_time;
    int K = params.Rw.k;

    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_simluation_by_lambda";
    Utils::createDirectory(params.OutputDir);
    std::string filename_mean_std_time = params.OutputDir + "/time_std_values.txt";

    std::ofstream arq;
    
    for (int i = 0; i < params.Rw.lambda_vector.size(); i++) {
        params.Rw.lambda = params.Rw.lambda_vector[i];
        // phase transition
        if (avg_epidemic_time >= params.Time) {
            std::cout << "Phase transition achieved, stop varying lambda" << std::endl;
            
            avg_epidemic_time = params.Time;
            std_epidemic_time = 0.012;
        }
        else
            start_simulation(params, strParams);
        
        lambda_mean_time.push_back(std::make_pair(params.Rw.lambda_vector[i], avg_epidemic_time));
        lambda_std_time.push_back(std::make_pair(params.Rw.lambda_vector[i], std_epidemic_time));
        
        arq.open(filename_mean_std_time, std::ofstream::app);
        arq << Utils::replace_dot_to_comma(params.Rw.lambda) << ";" << Utils::replace_dot_to_comma(avg_epidemic_time) << ";" << Utils::replace_dot_to_comma(std_epidemic_time) << "\n";
        arq.close();
    }

    EpidemicAnalysis * ep = new EpidemicAnalysis(params.OutputDir);
    ep->params = params;
    ep->mean_time_of_epidemic_over_n(lambda_mean_time, K);
    ep->std_time_of_epidemic_over_n(lambda_std_time, K);
    delete ep;
}

void EpidemicManager::start_simulation(Params params, std::string strParams) {
    epidemic_time = 0;
    std_epidemic_time = 0;
    means.clear();

    ManipulaGrafoV graph = create_graph(params);

    std::string fName = Utils::datetime_now_to_string();

    params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Runs);
    +"_" + std::to_string(params.Runs);
    output_dir = params.OutputDir;
    Utils::createDirectory(params.OutputDir);

    std::ofstream arq;
    std::string paramsFile = params.OutputDir + "/Params_Execution.json";

    arq.open(paramsFile);
    arq << strParams << "\n";
    arq.close();

    std::cout << "Starting simulation" << "\n";

    int num_threads = params.Runs;

    //std::thread * threads = new std::thread[num_threads];

#pragma omp parallel for
    for (int i = 0; i < params.Runs; i++) {
        //threads[i] = std::thread(&EpidemicManager::run_simulation_in_thread, this, params, strParams, graph, i);
        run_simulation_in_thread(params, strParams, graph, i);
    }

    //for (int i = 0; i < params.Runs; i++)
    //        threads[i].join();

    avg_epidemic_time = epidemic_time / (double) params.Runs;
    std_epidemic_time = std_desviation(avg_epidemic_time);

    std::cout << "\n";
    std::cout << "N: " << params.Graph.n << "\n";
    std::cout << "K: " << params.Rw.k << "\n";
    std::cout << "Tempo médio de execução da epidemia " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << std::endl;
    std::cout << "Tempo médio de simulação da epidemia " << avg_epidemic_time << std::endl;
    std::cout << "\n";

    std::string file_all = params.OutputDir + "/results.txt";

    arq.open(file_all, std::ofstream::out | std::ofstream::app);
    //arq << "Tempo médio da epidemia: " << boost::accumulators::mean(time_epidemic_acc) << '\n';
    arq << "Tempo médio da epidemia: " << avg_epidemic_time << '\n';
    //arq << "Desvio padrão tempo da epidemia: " << sqrt(boost::accumulators::variance(time_epidemic_acc)) << " | " << std_epidemic_time << '\n';
    arq << "Desvio padrão tempo da epidemia: " << std_epidemic_time << '\n';

    arq << "Tempo médio de execução: " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << '\n';
    // don't work when used with mean_time_epidemic... accumulator doesn't have clear method
    //arq << "Desvio padrão tempo de execução: " << sqrt(boost::accumulators::variance(time_epidemic_execution_acc)) / 1000. << "s" << '\n';

    arq.close();

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
    paramsR.OutputDir = params.OutputDir + "/Round " + std::to_string(run);
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
    time_epidemic_execution_acc(duration);
    means_dur_execution.push_back(duration);

    epidemic_time += sim->time;
    time_epidemic_acc(sim->time);
    means.push_back(sim->time);
    mtx.unlock();
    double time_s = duration / 1000.;

    std::cout << "Tempo de execução da epidemia: " << time_s << "s" << std::endl;

    std::ofstream arq;
    arq.open(sim->file_name_system_results, std::ofstream::out | std::ofstream::app);
    arq << "Tempo da epidemia: " << sim->time << '\n';
    arq << "Tempo de execução: " << time_s << "s" << '\n';
    arq.close();

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

double EpidemicManager::std_desviation(double mean) {
    double sum = 0.0;
    for (int i = 0; i < means.size(); i++)
        sum += std::pow(means.at(i) - mean, 2);
    sum = sum / (double) means.size();
    return sqrt(sum);
}