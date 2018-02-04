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
    Logger::Info("Running epidemic changing k");
    
    std::vector<std::pair<double, double> > k_mean_time;
    std::vector<std::pair<double, double> > k_std_time;
    int N = params.Graph.n;

    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_simluation_by_k_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.ki) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Rw.gama) + "_" + std::to_string(params.Runs);
    Utils::createDirectory(params.OutputDir);
    std::string filename_mean_std_time = params.OutputDir + "/time_std_values.txt";

    std::ofstream arq;

    for (int i = 0; i < params.KS_vector.size(); i++) {
        params.Rw.k = params.KS_vector[i];
        params.Rw.ki = params.Rw.kiPer * params.Rw.k < 1 ? 1 : params.Rw.kiPer * params.Rw.k;
        
        // phase transition
        if (avg_epidemic_time >= params.Time) {
            Logger::Info("Phase transition achieved, stop varying k. Setting Time limit for avg_time");

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
    Logger::Info("Running epidemic changing n");
    
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
    Logger::Info("Running epidemic changing lambda");
    
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
    Logger::Info("Running epidemic.\nParameters:\n   graph: " 
            + std::to_string(params.Graph.Type) 
            + "\n   n: " + std::to_string(params.Graph.n)
            + "\n   k: " + std::to_string(params.Rw.k)
            + "\n   ki: " + std::to_string(params.Rw.kiPer)
            + "\n   lambda: " + std::to_string(params.Rw.lambda)
            + "\n   gama: " + std::to_string(params.Rw.gama)
            + "\n   runs: " + std::to_string(params.Runs));
    
    epidemic_time = 0;
    std_epidemic_time = 0;
    means.clear();

    ManipulaGrafoV graph = create_graph(params);

    std::string fName = Utils::datetime_now_to_string();

    params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Runs);
    +"_" + std::to_string(params.Runs);
    output_dir = params.OutputDir;
    Utils::createDirectory(params.OutputDir);

    std::ofstream arq;
    std::string paramsFile = params.OutputDir + "/Params_Execution.json";

    arq.open(paramsFile);
    arq << strParams << "\n";
    arq.close();

    int num_threads = params.Runs;

    //std::thread * threads = new std::thread[num_threads];

#pragma omp parallel for
    for (int i = 0; i < params.Runs; i++) {
        run_simulation_in_thread(params, strParams, graph, i);
    }

    avg_epidemic_time = epidemic_time / (double) params.Runs;
    TimeStatistics time_stats = time_epidemic_statistics(avg_epidemic_time);
    std_epidemic_time = time_stats.std;

    /*std::cout << "\n";
    std::cout << "N: " << params.Graph.n << "\n";
    std::cout << "K: " << params.Rw.k << "\n";
    std::cout << "Tempo médio de execução da epidemia " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << std::endl;
    std::cout << "Tempo médio de simulação da epidemia " << avg_epidemic_time << std::endl;
    std::cout << "\n";*/

    std::string file_all = params.OutputDir + "/results.txt";

    arq.open(file_all, std::ofstream::out | std::ofstream::app);
    arq << "Parametros\n";
    arq << "K: " << params.Rw.k << "\n";
    arq << "Ki_Per: " << params.Rw.kiPer << " Ki: " << params.Rw.ki << "\n";
    arq << "Tempo médio da epidemia: " << avg_epidemic_time << '\n';
    arq << "Desvio padrão tempo da epidemia: " << time_stats.std << '\n';
    arq << "Menor tempo: " << time_stats.min << "\n";
    arq << "Maior tempo: " << time_stats.max << "\n";
    arq << "Tempo médio de execução: " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << '\n';
    arq << "\n";
    arq << "Rodada" << std::setw(20) << "Tempo" << "\n";
    for(int i=0; i < means.size(); i++)
        arq << std::get<0>(means.at(i)) << std::setw(25) << std::get<1>(means.at(i)) << "\n";

    arq.close();
}

void EpidemicManager::run_simulation_in_thread(Params params, std::string strParams, ManipulaGrafoV graph, int run) {
    Logger::Info("--------- ROUND " + std::to_string(run) + " --------- \n");

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
    means.push_back(std::make_tuple(run, sim->time));
    mtx.unlock();
    double time_s = duration / 1000.;

    //std::cout << "Tempo de execução da epidemia: " << time_s << "s" << std::endl;

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
            ep->getRandomWalkCCDFAndStatistics(sim->randomWalks.at(i)->file_name_sample_results, std::to_string(i), sim->randomWalks.at(i)->gama, sim->randomWalks.at(i)->tau);
            ep->RandomWalkWalkingCCDF(sim->randomWalks.at(i)->file_name_walking_times, std::to_string(i), sim->randomWalks.at(i)->lambda);
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

    Logger::Info("--------- END ROUND " + std::to_string(run) + " --------- \n");

}

TimeStatistics EpidemicManager::time_epidemic_statistics(double mean) {
    double sum = 0.0;
    double max = 0.0;
    double min = INT64_MAX;
    for (int i = 0; i < means.size(); i++)
    {
        double time = std::get<1>(means.at(i));
        sum += std::pow(time - mean, 2);
        if(time > max)
            max = time;
        if(time < min)
            min = time;
    }
    
    sum = sum / (double) means.size();
    double std = sqrt(sum);
    TimeStatistics t;
    t.mean = mean;
    t.max = max;
    t.min = min;
    t.std = std;
    
    return t;
}