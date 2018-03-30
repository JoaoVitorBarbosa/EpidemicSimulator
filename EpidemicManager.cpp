/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "EpidemicManager.h"

EpidemicManager::EpidemicManager() {
    aggregate_time_ms = 0;
    epidemic_time = 0.0;
    do_analysis = true;
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
    Logger::Info("mean_time_epidemic_by_k");

    std::vector<std::pair<double, double> > k_mean_time;
    std::vector<std::pair<double, double> > k_std_time;
    int N = params.Graph.n;

    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_simulation_by_k_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.ki) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Rw.gama) + "_" + std::to_string(params.Runs);
    Utils::createDirectory(params.OutputDir);
    std::string filename_mean_std_time = params.OutputDir + "/time_std_values.txt";

    std::ofstream arq;

    for (int i = 0; i < params.KS_vector.size(); i++) {
        params.Rw.k = params.KS_vector[i];
        params.Rw.ki = params.Rw.kiPer * params.Rw.k < 1 ? 1 : params.Rw.kiPer * params.Rw.k;

        // phase transition
        if (avg_epidemic_time >= params.Time) {
            Logger::Info("K= " + std::to_string(params.Rw.k) + " Phase transition achieved.. Setting Time limit for avg_time");

            avg_epidemic_time = params.Time;
            std_epidemic_time = 0.0012;
        } else {
            // don't use parallelism
            start_simulation(params, strParams, true);
        }

        //k_mean_time.push_back(std::make_pair(params.KS_vector[i], avg_epidemic_time));
        //k_std_time.push_back(std::make_pair(params.KS_vector[i], std_epidemic_time));

        arq.open(filename_mean_std_time, std::ofstream::app);
        arq << params.Rw.k << "," << avg_epidemic_time << "," << std_epidemic_time << "\n";
        arq.close();
    }

    if (do_analysis) {
        EpidemicAnalysis * ep = new EpidemicAnalysis(params.OutputDir);
        ep->params = params;
        ep->mean_time_of_epidemic_over_k(k_mean_time, N);
        ep->std_time_of_epidemic_over_k(k_std_time, N);
        delete ep;
    }

}

void EpidemicManager::mean_time_epidemic_by_n(Params params, std::string strParams) {
    Logger::Info("Running epidemic changing n");

    std::vector<std::pair<double, double> > n_mean_time;
    std::vector<std::pair<double, double> > n_std_time;
    int K = params.Rw.k;
    std::string outputdir;
    for (int i = 0; i < params.Graph.NS_vector.size(); i++) {
        params.Graph.n = params.Graph.NS_vector[i];
        start_simulation(params, strParams, true);
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
            Logger::Info("Phase transition achieved, stop varying lambda");

            avg_epidemic_time = params.Time;
            std_epidemic_time = 0.012;
        } else
            start_simulation(params, strParams, true);

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

void EpidemicManager::measure_scalability_by_random_walks(Params params, std::string strParams) {
    Logger::Info("Scalability by random walks");

    std::vector<std::pair<double, double> > mean_time;
    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_scalability_by_k_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.ki) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Rw.gama) + "_" + std::to_string(params.Runs);
    Utils::createDirectory(params.OutputDir);
    std::string mean_time_file = params.OutputDir + "/scalability_by_k.txt";

    std::ofstream arq;

    for (int i = 0; i < params.KS_vector.size(); i++) {
        params.Rw.k = params.KS_vector[i];
        params.Rw.ki = params.Rw.kiPer * params.Rw.k < 1 ? 1 : params.Rw.kiPer * params.Rw.k;

        start_simulation(params, strParams, false);

        Logger::Info("Tempo de Epidemia: " + std::to_string(avg_epidemic_time));

        double mean_duration = aggregate_time_ms / params.Runs;
        double m_duration = mean_duration / 1000.0;
        mean_time.push_back(std::make_pair(params.Rw.k, m_duration));
        TimeStatistics stats = duration_standard_desviation(means_dur_execution, mean_duration);

        arq.open(mean_time_file, std::ofstream::app);
        arq << params.Rw.k << "," << m_duration << "," << stats.std / 1000.0 << "," << avg_epidemic_time << "\n";
        arq.close();
    }
}

void EpidemicManager::measure_scalability_by_vertices(Params params, std::string strParams) {
    Logger::Info("Scalability by vertices");

    std::vector<std::pair<double, double> > mean_time;
    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_scalability_by_vertices_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.ki) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Rw.gama) + "_" + std::to_string(params.Runs);
    Utils::createDirectory(params.OutputDir);
    std::string mean_time_file = params.OutputDir + "/scalability_by_vertices.txt";

    std::ofstream arq;

    for (int i = 0; i < params.Graph.NS_vector.size(); i++) {
        params.Graph.n = params.Graph.NS_vector[i];
        params.Rw.ki = params.Rw.kiPer * params.Rw.k < 1 ? 1 : params.Rw.kiPer * params.Rw.k;

        start_simulation(params, strParams, false);

        Logger::Info("Tempo de Epidemia: " + std::to_string(avg_epidemic_time));

        double mean_duration = aggregate_time_ms / params.Runs;
        double m_duration = mean_duration / 1000.0;
        mean_time.push_back(std::make_pair(params.Graph.n, m_duration));
        TimeStatistics stats = duration_standard_desviation(means_dur_execution, mean_duration);

        arq.open(mean_time_file, std::ofstream::app);
        arq << params.Graph.n << "," << m_duration << "," << stats.std / 1000.0 << "," << avg_epidemic_time << "\n";
        arq.close();
    }
}

void EpidemicManager::measure_scalability_by_time(Params params, std::string strParams) {
    Logger::Info("Scalability by time");

    std::vector<std::pair<double, double> > mean_time;
    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_scalability_by_time_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.ki) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Rw.gama) + "_" + std::to_string(params.Runs);
    Utils::createDirectory(params.OutputDir);
    std::string mean_time_file = params.OutputDir + "/scalability_by_time.txt";

    std::ofstream arq;

    for (int i = 0; i < params.Time_vector.size(); i++) {
        params.Time = params.Time_vector[i];
        params.Rw.ki = params.Rw.kiPer * params.Rw.k < 1 ? 1 : params.Rw.kiPer * params.Rw.k;

        start_simulation(params, strParams, false);

        Logger::Info("Tempo de Epidemia: " + std::to_string(avg_epidemic_time));

        double mean_duration = aggregate_time_ms / params.Runs;
        double m_duration = mean_duration / 1000.0;
        mean_time.push_back(std::make_pair(params.Time, m_duration));
        TimeStatistics stats = duration_standard_desviation(means_dur_execution, mean_duration);

        arq.open(mean_time_file, std::ofstream::app);
        arq << params.Time << "," << m_duration << "," << stats.std / 1000.0 << "," << avg_epidemic_time << "\n";
        arq.close();
    }
}

void EpidemicManager::start_simulation(Params params, std::string strParams, bool do_parallel) {
    Logger::Info("Running epidemic.\nParameters:\n   graph: "
            + std::to_string(params.Graph.Type)
            + "\n   Time: " + std::to_string(params.Time)
            + "\n   n: " + std::to_string(params.Graph.n)
            + "\n   k: " + std::to_string(params.Rw.k)
            + "\n   ki: " + std::to_string(params.Rw.ki)
            + "\n   lambda: " + std::to_string(params.Rw.lambda)
            + "\n   gama: " + std::to_string(params.Rw.gama)
            + "\n   runs: " + std::to_string(params.Runs));

    epidemic_time = 0;
    std_epidemic_time = 0;
    aggregate_time_ms = 0;
    means.clear();
    means_dur_execution.clear();

    ManipulaGrafoV graph = create_graph(params);

    std::string fName = Utils::datetime_now_to_string();
    params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + params.OutputDir = params.OutputDir + "/" + fName + "_" + std::to_string(params.Time) + "_" + std::to_string(params.Graph.n) + "_" + std::to_string(params.Rw.k) + "_" + std::to_string(params.Rw.kiPer) + "_" + std::to_string(params.Rw.lambda) + "_" + std::to_string(params.Runs);
    +"_" + std::to_string(params.Runs);
    output_dir = params.OutputDir;
    Utils::createDirectory(params.OutputDir);

    Logger::Trace("Output Directory created: " + output_dir);

    std::ofstream arq;
    std::string paramsFile = params.OutputDir + "/Params_Execution.json";

    arq.open(paramsFile);
    arq << strParams << "\n";
    arq.close();

    if (do_parallel) {
        Logger::Trace("Starting Simulation parallel");
#pragma omp parallel for
        for (int i = 0; i < params.Runs; i++) {
            run_simulation_in_thread(params, strParams, graph, i);
        }
    } else {
        Logger::Trace("Starting Simulation serial");
        for (int i = 0; i < params.Runs; i++) {
            run_simulation_in_thread(params, strParams, graph, i);
        }
    }

    avg_epidemic_time = epidemic_time / (double) params.Runs;
    TimeStatistics time_stats = time_epidemic_statistics(avg_epidemic_time);
    std_epidemic_time = time_stats.std;

    std::string file_all = params.OutputDir + "/results.txt";

    arq.open(file_all, std::ofstream::out | std::ofstream::app);
    arq << "Parametros\n";
    arq << "N: " << params.Graph.n << "\n";
    arq << "Grafo: " << params.Graph.Type << "\n";
    arq << "K: " << params.Rw.k << "\n";
    arq << "Ki_Per: " << params.Rw.kiPer << " Ki: " << params.Rw.ki << "\n";
    arq << "Tempo médio da epidemia: " << avg_epidemic_time << '\n';
    arq << "Desvio padrão tempo da epidemia: " << time_stats.std << '\n';
    arq << "Menor tempo: " << time_stats.min << "\n";
    arq << "Maior tempo: " << time_stats.max << "\n";
    arq << "Tempo médio de execução: " << (aggregate_time_ms / params.Runs) / 1000.00 << "s" << '\n';
    arq << "\n";
    arq << "Rodada" << std::setw(20) << "Tempo" << "\n";
    for (int i = 0; i < means.size(); i++)
        arq << std::get<0>(means.at(i)) << std::setw(25) << std::get<1>(means.at(i)) << "\n";

    arq.close();

    //if (do_analysis) {
        Logger::Trace("Infected avg over time");
        std::string analysisDir = params.OutputDir;
        EpidemicAnalysis * ep = new EpidemicAnalysis(analysisDir);
        ep->params = params;
        ep->infected_in_time_in_runs(params.Runs, params.OutputDir);

        delete ep;
    //}
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
    // get final time
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    // calculate duration of epidemic
    int duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();

    begin_time = std::chrono::high_resolution_clock::now();
    sim->process();
    end_time = std::chrono::high_resolution_clock::now();
    int duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();

    Logger::Trace("Time building simulator: " + std::to_string(duration1 / 1000.));
    Logger::Trace("Time processing events: " + std::to_string(duration / 1000.));
    Logger::Trace("Tempo de epidemia: " + std::to_string(sim->time));

    std::mutex mtx;
    mtx.lock();
    aggregate_time_ms += duration;
    time_epidemic_execution_acc(duration);
    means_dur_execution.push_back(duration);
    epidemic_time += sim->time;
    means.push_back(std::make_tuple(run, sim->time));
    mtx.unlock();
    double time_s = duration / 1000.;

    //std::cout << "Tempo de execução da epidemia: " << time_s << "s" << std::endl;

    std::ofstream arq;
    arq.open(sim->file_name_system_results, std::ofstream::out | std::ofstream::app);
    arq << "Tempo de epidemia: " << sim->time << '\n';
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
//            ep->stateDistribution(sim->vertices.at(i)->fileNameTimeResult,  std::to_string(i));

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
    for (int i = 0; i < means.size(); i++) {
        double time = std::get<1>(means.at(i));
        sum += std::pow(time - mean, 2);
        if (time > max)
            max = time;
        if (time < min)
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

TimeStatistics EpidemicManager::duration_standard_desviation(std::vector<double> times, double mean) {
    double sum = 0.0;
    double max = 0.0;
    double min = INT64_MAX;
    for (int i = 0; i < times.size(); i++) {
        double time = times.at(i);
        sum += std::pow(time - mean, 2);
        if (time > max)
            max = time;
        if (time < min)
            min = time;
    }

    sum = sum / (double) times.size();
    double std = sqrt(sum);
    TimeStatistics t;
    t.mean = mean;
    t.max = max;
    t.min = min;
    t.std = std;

    return t;
}