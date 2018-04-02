/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "EpidemicAnalysis.h"
#include <set>
#include <math.h>  

EpidemicAnalysis::EpidemicAnalysis(std::string _outputDir) {
    outputDir = _outputDir;
    do_validation = true;
}

/*
 * 
    Privates
 * 
 */

std::vector<std::pair<double, double> > EpidemicAnalysis::getCCDFPoints(std::multiset<double> elements) {
    std::vector<std::pair<double, double> > xy_pts;
    int i = 0;
    int size = elements.size();
    for (std::multiset<double>::iterator it = elements.begin(); it != elements.end(); ++it) {
        double frac = (size - i) / (double) size;
        xy_pts.push_back(std::make_pair(*it, frac));
        i++;
    }

    return xy_pts;
}

void EpidemicAnalysis::createStateStatisticsAndCCDF(std::vector<std::pair<double, int> > timeStampStatesChange, std::string rw, double gama, double tau) {
    std::multiset<double> suscetibles_segments;
    std::multiset<double> contracteds_segments;
    std::multiset<double> infecteds_segments;

    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::variance>> intervalsInfected;
    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::variance>> intervalsContracted;
    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::variance>> intervalsSusceptible;

    // states: 0 -> Susceptible, 1 -> Contracted, 2 -> infected
    // begins from 1 because first element is always (zero, initial state)
    for (int i = 1; i < timeStampStatesChange.size(); i++) {
        // avoid to create negative intervals when parsing all random walks
        if (timeStampStatesChange[i].first == 0)
            continue;

        double interval = timeStampStatesChange[i].first - timeStampStatesChange[i - 1].first;

        if (interval < 0)
            std::cout << interval << std::endl;

        // time infected
        if (timeStampStatesChange[i].second == 0) {
            //std::cout << interval << std::endl;
            infecteds_segments.insert(interval);
            intervalsInfected(interval);
        }
        // time suscetible
        if (timeStampStatesChange[i].second == 1) {
            suscetibles_segments.insert(interval);
            intervalsSusceptible(interval);
        }
        // time contracted
        if (timeStampStatesChange[i].second == 2) {
            contracteds_segments.insert(interval);
            intervalsContracted(interval);
        }
    }

    std::vector<std::pair<double, double> > susceptiblesPts = getCCDFPoints(suscetibles_segments);
    std::vector<std::pair<double, double> > contractedsPts = getCCDFPoints(contracteds_segments);
    std::vector<std::pair<double, double> > infectedsPts = getCCDFPoints(infecteds_segments);

    std::string prefix(rw + "_");
    std::string filename = this->outputDir + "/" + prefix + "_CCDFEstadoEpidemico.png";

    // print all ccdf in one graphic
    Gnuplot gp;

    if (do_validation) {
        //gp << "set title 'CCDF tempo gasto no estado infectado' \n";
        gp << "set encoding utf8 \n";
        gp << "set xlabel 'Tempo de permanência no estado Infectado' font 'Helvetica,18' \n";
        gp << "set ylabel 'P[T>=t]' font 'Helvetica,18' \n";
        gp << "set term png \n";
        gp << "set output '" + this->outputDir + "/" + prefix + "CCDFInfectado.png" + "' \n";
        gp << "set grid \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(infectedsPts) << "title 'Amostral' with lines lt rgb '#00FF00', exp(-" + std::to_string(gama) + "*x) title 'Estatístico' with lines lt rgb '#FF0000' \n";

        //gp << "set title 'CCDF tempo gasto no estado contraido' \n";
        gp << "set encoding utf8 \n";
        gp << "set xlabel 'Tempo de permanência  no estado contraído' font 'Helvetica,18' \n";
        gp << "set ylabel 'P[T>=t]' font 'Helvetica,18' \n";
        gp << "set term png \n";
        gp << "set output '" + this->outputDir + "/" + prefix + "_CCDFContraido.png" + "' \n";
        gp << "set grid \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(contractedsPts) << "title 'Amostral' with lines lt rgb '#00FF00', exp(-" + std::to_string(tau) + "*x) title 'Estatístico' with lines lt rgb '#FF0000' \n";

        //gp << "set title 'CCDF tempo gasto no estado \n";
        gp << "set encoding utf8 \n";
        gp << "set xlabel 'Tempo de permanência  no estado epidêmico' font 'Helvetica,18' \n";
        gp << "set ylabel 'P[T>=t]' font 'Helvetica,18' \n";
        gp << "set term png \n";
        gp << "set output '" + filename + "' \n";
        gp << "set grid \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(susceptiblesPts) << "title 'Suscetível' with lines lt rgb '#00FF00', "
                << gp.file1d(contractedsPts) << "title 'Contraído' with lines lt rgb '#0000FF', "
                << gp.file1d(infectedsPts) << "title 'Infectado' with lines lt rgb '#FF0000' \n";
    }

    std::ofstream arq;
    arq.open(this->outputDir + "/" + prefix + "_Results.txt", std::ofstream::out | std::ofstream::app);

    double infected_mean = boost::accumulators::mean(intervalsInfected);
    double sum = 0;
    for(int i = 0; i < infectedsPts.size(); i++)
        sum += infectedsPts[i].second;

    arq << "Estado Infectado \n"
            << "    Média teórica:          " << 1 / gama << "    Média empírica: " << infected_mean << "\n"
            << "    Desvio Padrão teórico:  " << 1 / gama << "    Desvio padrão empírico: " << std::sqrt(boost::accumulators::variance(intervalsInfected)) << "\n"
            << "    Mediana teórica:        " << std::log(2) / gama << "    Mediana empírica: " << get_median(infecteds_segments).median << "\n"
            << "    Gama teórico:           " << gama << "    Gama empírico: " << 1 / infected_mean << "\n"
            << "    Gama teórico:           " << gama << "    Gama empírico-2: " << infectedsPts.size() / sum << "\n"
            << '\n';

    arq << "Estado Contraído \n"
            << "    Média teórica:          " << 1 / tau << "    Média empírica: " << boost::accumulators::mean(intervalsContracted) << "\n"
            << "    Desvio Padrão teórico:  " << 1 / tau << "    Desvio padrão empírico: " << std::sqrt(boost::accumulators::variance(intervalsContracted)) << "\n"
            << "    Mediana teórica:        " << std::log(2) / tau << "    Mediana empírica: " << get_median(contracteds_segments).median << "\n"
            << "    Tau teórico:           " << tau << "    Tau empírico: " << 1 / boost::accumulators::mean(intervalsContracted) << "\n"
            << '\n';

    arq << "Estado Suscetível \n"
            << "    Média: " << boost::accumulators::mean(intervalsSusceptible) << "\n"
            << "    Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsSusceptible)) << "\n"
            << "    Mediana: " << get_median(suscetibles_segments).median << "\n"
            << '\n';

    arq.close();
}

std::vector<std::pair<int, double>> EpidemicAnalysis::getStateInterval(std::string filepath) {
    std::vector<std::pair<int, double>> timeStampStatesChange;

    std::ifstream arq;
    arq.open(filepath.c_str());

    std::string linha = "";

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        if (line.find('#') == 0)
            continue;

        int p = line.find_first_of(',');
        if (p == -1)
            break;

        int rw_inf = atof(line.substr(0, p).c_str());
        int next = p + 1;
        std::string str = line.substr(next, line.length());
        double interval = atof(str.c_str());

        std::pair<int, double> pair(rw_inf, interval);
        timeStampStatesChange.push_back(pair);
    }

    arq.close();

    return timeStampStatesChange;
}


/* ---------------------------------------------------------------------------------------------- */

/*
    Publics
 */

void EpidemicAnalysis::readTimestampStateChangeCSV(std::string filename, std::vector<std::pair<double, int> > &timeStampStatesChange) {
    std::ifstream arq;
    arq.open(filename.c_str());

    std::string linha = "";

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        if (line.find('#') == 0)
            continue;

        int p = line.find_first_of(',');
        if (p == -1)
            break;

        double time = atof(line.substr(0, p).c_str());
        int next = p + 1;
        std::string str = line.substr(next, line.length());
        int state = atof(str.c_str());

        std::pair<double, int> pair(time, state);
        timeStampStatesChange.push_back(pair);
    }

    arq.close();
}

void EpidemicAnalysis::getSystemCCDFAndStatistics(std::string dir, int num_rw) {
    std::vector<std::pair<double, int> > all;
    for (int i = 0; i < num_rw; i++) {
        readTimestampStateChangeCSV(dir + "/RW" + std::to_string(i) + "_Results.txt", all);
        all.push_back(std::make_pair(0, 0));
    }

    //createStateStatisticsAndCCDF(all, "System");
}

void EpidemicAnalysis::readWalkingTime(std::string filename, std::multiset<double> &timeStampWalking) {
    std::ifstream arq;
    arq.open(filename.c_str());

    std::string linha = "";
    double lastTime = 0.0;

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        if (line.find('#') == 0)
            continue;

        double time = atof(line.c_str());
        double interval = time - lastTime;
        if (interval > 0)
            timeStampWalking.insert(interval);

        lastTime = time;
    }

    arq.close();
}

void EpidemicAnalysis::RandomWalkWalkingCCDF(std::string filepath, std::string randomWalk, double lambda) {
    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::variance>> intervalsWalking;
    std::multiset<double> walking_segments;
    readWalkingTime(filepath, walking_segments);

    intervalsWalking = std::for_each(walking_segments.begin(), walking_segments.end(), intervalsWalking);

    std::vector<std::pair<double, double> > walkingPts = getCCDFPoints(walking_segments);

    std::string filename = this->outputDir + "/" + randomWalk + "_CCDF_caminhar.png";

    if (do_validation) {

        // writing in file
        std::string prefix(randomWalk + "_");

        double emp_mean = boost::accumulators::mean(intervalsWalking);
        Stats stats = get_median(walking_segments);
        
        std::ofstream arq;
        arq.open(this->outputDir + "/" + prefix + "_Results.txt", std::ofstream::out | std::ofstream::app);

        arq << "Dados do evento caminhar \n"
                << "    Média esperada:          " << 1 / lambda << "         Média empírico: " << emp_mean << "\n"
                << "    Desvio Padrão esperado:  " << 1 / lambda << "         Desvio padrão empírico: " << std::sqrt(boost::accumulators::variance(intervalsWalking)) << "\n"
                << "    Mediana teórica:        " << std::log(2) / lambda << "         Mediana empírico: " << stats.median << "\n"
                << "    Lambda teórico:           " << lambda << "         Gama empírico: " << 1 / emp_mean << "\n"
                << '\n';

        Gnuplot gp;

        //gp << "set title 'CCDF do tempo de caminhada\n";
        gp << "set encoding utf8 \n";
        gp << "set xlabel 'Tempo de permanência no vértice' font 'Helvetica,16' \n";
        gp << "set ylabel 'P[t>=T]' font 'Helvetica,18' \n";
        gp << "set term png \n";
        gp << "set grid \n";
        gp << "set output '" + filename + "' \n";
        gp << "set logscale y \n";

        gp << "plot" << gp.file1d(walkingPts) << "title 'Amostral' with lines lt rgb '#00FF00', exp(-" + std::to_string(lambda) + "*x) title 'Estatístico' with lines lt rgb '#FF0000' \n";

        //gp << "plot" << gp.file1d(walkingPts) << "title 'empiric' with lines lt rgb '#00FF00 \n'";
    }

}

void EpidemicAnalysis::randomWalkStateTimeSeries(std::string arquivo, int k, int time_total) {

    Logger::Trace("Gerando fração de indivíduos em cada estado");
    
    std::vector<std::pair<double, int> > infecteds;
    std::vector<std::pair<double, int> > contracteds;
    std::vector<std::pair<double, int> > susceptibles;

    std::string file_frac_infec = this->outputDir + "/frac_infected.csv";
    std::string file_frac_cont = this->outputDir + "/frac_contracted.csv";
    std::string file_frac_sus = this->outputDir + "/frac_susceptible.csv";
    
    std::ofstream arq_inf, arq_cont, arq_sus;
    arq_inf.open(file_frac_infec, std::ofstream::out | std::ofstream::app);
    arq_cont.open(file_frac_cont, std::ofstream::out | std::ofstream::app);
    arq_sus.open(file_frac_sus, std::ofstream::out | std::ofstream::app);
    
    std::ifstream arq;
    arq.open(arquivo.c_str());

    std::string linha = "";
    int max = 0;

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        if (line.find('#') == 0)
            continue;

        std::stringstream lineStream(line);
        std::string cell;

        int i = 0;
        double time;
        int inf, sus, cont;
        while (std::getline(lineStream, cell, ',')) {
            if (i == 0)
                time = atof(cell.c_str());
            if (i == 1) {
                inf = atof(cell.c_str());
                max = inf > max ? inf : max;
            }
            if (i == 2) {
                cont = atof(cell.c_str());
                max = cont > max ? cont : max;
            }
            if (i == 3) {
                sus = atof(cell.c_str());
                max = sus > max ? sus : max;
            }
            i++;
        }

//        Logger::Trace("Time: " + std::to_string(time));
//        Logger::Trace("Fração inf: " + std::to_string(inf/(double)k));
//        Logger::Trace("Fração cont: " + std::to_string(cont/(double)k));
//        Logger::Trace("Fração sus: " + std::to_string(sus/(double)k));
        
        arq_inf << time << ",      " << inf/(double)k << std::endl;
        arq_cont << time << ",      " << cont/(double)k << std::endl;
        arq_sus << time << ",      " << sus/(double)k << std::endl;
        
        infecteds.push_back(std::make_pair(time, inf));
        contracteds.push_back(std::make_pair(time, cont));
        susceptibles.push_back(std::make_pair(time, sus));
    }

    arq_inf.close();
    arq_cont.close();
    arq_sus.close();
    arq.close();

    std::string filename = this->outputDir + "/rwStateTimeSeries.png";

    Gnuplot gp;
    gp << "set title '' \n";
    gp << "set xlabel 't' \n";
    gp << "set ylabel 'Fracaoo de individuos' \n";
    gp << "set xrange [0:" + std::to_string(time_total) + "]\n";
    gp << "set yrange [0:" + std::to_string(k) + "]\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set terminal png medium size 1400,900 \n";
    gp << "plot '-' title 'Suscetivel' with lines lt rgb '#00FF00', '-' title 'Contraido' with lines lt rgb '#0000FF', '-' title 'Infectado' with lines lt rgb '#FF0000' \n";

    gp.send1d(susceptibles);
    gp.send1d(contracteds);
    gp.send1d(infecteds);

    Logger::Trace("End");
}

void EpidemicAnalysis::infected_in_time_in_runs(int runs, std::string output, int time, int k) {
    Logger::Trace("Begin: Fração média de indivíduos em cada estado");
    
    std::string file_path = outputDir + "/fracao_estado_no_tempo.csv";
    std::vector<std::pair<double, double> > infecteds_1;
    std::vector<std::pair<double, double> > infecteds_avg;
    std::vector<std::pair<double, double> > contracted_avg;
    std::vector<std::pair<double, double> > susceptible_avg;

    for (int round = 0; round < runs; round++) {
        std::string file_input = output + "/Round " + std::to_string(round) + "/system_results.txt";
        std::vector<std::pair<double, double> > infecteds_aux;
        std::vector<std::pair<double, double> > contracted_aux;
        std::vector<std::pair<double, double> > susceptible_aux;
        
        // gets data
        std::ifstream arq;
        arq.open(file_input.c_str());

        std::string linha = "";
        int max = 0;
        int count_to_avg = 1;
        int sum_inf = 0;
        int sum_cont = 0;
        int sum_sus = 0;
        while (!arq.eof()) {

            std::string line;
            std::getline(arq, line);

            if (line.find('#') == 0)
                continue;

            std::stringstream lineStream(line);
            std::string cell;

            int i = 0;
            double time;
            int inf, sus, cont;

            while (std::getline(lineStream, cell, ',')) {
                if (i == 0)
                    time = atof(cell.c_str());
                if (i == 1) {
                    inf = atof(cell.c_str());
                    max = inf > max ? inf : max;
                }
                if (i == 2) {
                    cont = atof(cell.c_str());
                    max = cont > max ? cont : max;
                }
                if (i == 3) {
                    sus = atof(cell.c_str());
                    max = sus > max ? sus : max;
                }
                i++;
            }

            sum_inf += inf;
            sum_cont += cont;
            sum_sus += sus;
            
            count_to_avg++;
            int d = 50; // slote de tempo para fazer a media
            if (count_to_avg % d == 0) {
                infecteds_aux.push_back(std::make_pair(time, sum_inf / (double) d));
                contracted_aux.push_back(std::make_pair(time, sum_cont / (double) d));
                susceptible_aux.push_back(std::make_pair(time, sum_sus / (double) d));

                int index = (count_to_avg / d) - 1;
                if (infecteds_avg.size() < index + 1)
                {
                    infecteds_avg.push_back(std::make_pair(time, sum_inf / (double) d));
                    contracted_avg.push_back(std::make_pair(time, sum_cont / (double) d));
                    susceptible_avg.push_back(std::make_pair(time, sum_sus / (double) d));
                }
                else {
                    std::pair<double, double> aux = infecteds_avg[index];
                    std::pair<double, double> aux2 = contracted_avg[index];
                    std::pair<double, double> aux3 = susceptible_avg[index];
                    if (aux.first == 0.0)
                        aux.first = time;
                    infecteds_avg[index] = std::make_pair(aux.first, aux.second + sum_inf / (double) d);
                    contracted_avg[index] = std::make_pair(aux2.first, aux2.second + sum_cont / (double) d);
                    susceptible_avg[index] = std::make_pair(aux3.first, aux3.second + sum_sus / (double) d);
                }

                sum_inf = 0;
                sum_cont = 0;
                sum_sus = 0;
            }
        }

        if (infecteds_aux.size() > infecteds_1.size())
            infecteds_1 = infecteds_aux;

        arq.close();
    }

    for (int i = 0; i < infecteds_avg.size(); i++) {
        std::pair<double, double> aux = infecteds_avg[i];
        infecteds_avg[i] = std::make_pair(aux.first, (aux.second / (double) runs)/ (double) k); // (/runs) média das rodadas 
        if(i < contracted_avg.size())
        {
            aux = contracted_avg[i];
            contracted_avg[i] = std::make_pair(aux.first, (aux.second / (double) runs)/ (double) k); // (/runs) média das rodadas 
        }
        
        if(i < susceptible_avg.size())
        {
            aux = susceptible_avg[i];
            susceptible_avg[i] = std::make_pair(aux.first, (aux.second / (double) runs)/ (double) k); // (/runs) média das rodadas 
        }
    }

    // remove ruido
    // infecteds_avg.pop_back();


    std::ofstream arq;
    arq.open(file_path, std::ofstream::out | std::ofstream::app);
    arq << "#time,uma rodada infectado,infectado,contradio,suscetivel" << std::endl;

    for (int i = 0; i < infecteds_avg.size(); i++) {
        arq << infecteds_avg[i].first << ",      " << infecteds_1[i].second << ",      " << infecteds_avg[i].second;
        if(i < contracted_avg.size())
        {
            arq << ",      " << contracted_avg[i].second;
        }
        
        if(i < susceptible_avg.size())
        {
            arq << ",      " << susceptible_avg[i].second;
        }
        
        arq << std::endl;
    }

    arq.close();

    std::string filename = this->outputDir + "/fracao_medio_estados.png";
    Gnuplot gp;
    gp << "set encoding utf8 \n";
    gp << "set title '' \n";
    gp << "set xlabel 'Tempo' \n";
    gp << "set ylabel 'I(t)' \n";
    gp << "set grid \n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set terminal png medium size 1400,900 \n";
    gp << "set yrange [0:1]\n";
    gp << "set xrange [0:" + std::to_string(time) + "]\n";
    gp << "plot '-' title 'Infectado' with lines lw 2 lt rgb '#FF0000', '-' title 'Contraido' with lines lw 2 lt rgb '#00FF00', '-' title 'Suscetivel' with lines lw 2 lt rgb '#0000FF' \n";
    //gp << "plot '-' title 'avg' with lines lw 10 lt rgb '#00FF00' \n";

    //gp.send1d(infecteds_1);
    gp.send1d(infecteds_avg);
    gp.send1d(contracted_avg);
    gp.send1d(susceptible_avg);
    
    Logger::Trace("End");
}

void readTimeNumberInfect(std::string filename, std::set<double, int>) {
    std::ifstream arq;
    arq.open(filename.c_str());

    std::string linha = "";

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        if (line.find('#') == 0)
            continue;

        int p = line.find_first_of(',');
        if (p == -1)
            break;

        double time = atof(line.substr(0, p).c_str());
        int next = p + 1;
        std::string str = line.substr(next, line.length());
        int number = atof(str.c_str());

        std::pair<double, int> pair(time, number);
        //timeStampStatesChange.push_back(pair);
    }

    arq.close();
}

void EpidemicAnalysis::getRandomWalkCCDFAndStatistics(std::string filepath, std::string title, double gama, double tau) {
    std::vector<std::pair<double, int> > pairs;
    readTimestampStateChangeCSV(filepath, pairs);
    createStateStatisticsAndCCDF(pairs, title, gama, tau);
}

void EpidemicAnalysis::stateDistribution(std::string filepathInfected, std::string filepathContracted, std::string filepathSusceptible, std::string title) {

    auto timeStampInfectedChange = getStateInterval(filepathInfected);
    auto timeStampContractedChange = getStateInterval(filepathContracted);
    auto timeStampSusceptibleChange = getStateInterval(filepathSusceptible);

    std::string filename = outputDir + "/StatesDistribution" + title + ".png";

    Gnuplot gp;

    gp << "set encoding utf8 \n";
    gp << "set xlabel 'Número de indivíduos' font 'Helvetica,16'\n";
    gp << "set ylabel 'Tempo(%)' font 'Helvetica,16'\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set grid \n";
    gp << "plot '-' title 'Suscetível' with lines lt rgb '#00FF00', '-' title 'Contraído' with lines lt rgb '#0000FF', '-' title 'Infectado' with lines lt rgb '#FF0000' \n";

    gp.send1d(timeStampSusceptibleChange);
    gp.send1d(timeStampContractedChange);
    gp.send1d(timeStampInfectedChange);

}

Stats EpidemicAnalysis::get_median(std::multiset<double> elems) {
    Stats stats;
    double median = 0.0;
    double mean = 0.0;
    int i = 0;
    int count = elems.size();
    for (auto it = elems.begin(); it != elems.end(); it++) {
        if (i == count / 2)
            median = *it;
        mean += *it;
        i++;
    }
    mean = mean / count;

    stats.mean = mean;
    stats.median = median;

    return stats;
}

void EpidemicAnalysis::mean_time_of_epidemic_over_k(std::vector<std::pair<double, double> > k_mean_time, int n) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_mean_time_by_k.gif";

    gp << "set title 'Tempo médio da epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Tempo médio' \n";
    gp << "set terminal gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'N=" + std::to_string(n) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_mean_time) << "title 'media' with linespoints lt rgb '#FF0000' \n";
}

void EpidemicAnalysis::std_time_of_epidemic_over_k(std::vector<std::pair<double, double> > k_std_time, int n) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_std_time_by_k.gif";

    gp << "set title 'Desvio padrão do tempo de epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Desvio padrão' \n";
    gp << "set termina gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'N=" + std::to_string(n) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_std_time) << "title 'desvio padrao' with linespoints lt rgb '#FF0000' \n";
}

void EpidemicAnalysis::mean_time_of_epidemic_over_n(std::vector<std::pair<double, double> > k_mean_time, int k) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_mean_time_by_n.gif";

    gp << "set title 'Tempo médio da epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Tempo médio' \n";
    gp << "set terminal gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'K=" + std::to_string(k) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_mean_time) << "title 'media' with linespoints lt rgb '#FF0000' \n";
}

void EpidemicAnalysis::std_time_of_epidemic_over_n(std::vector<std::pair<double, double> > k_std_time, int k) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_std_time_by_n.gif";

    gp << "set title 'Desvio padrão do tempo de epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Desvio padrão' \n";
    gp << "set termina gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'K=" + std::to_string(k) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_std_time) << "title 'desvio padrao' with linespoints lt rgb '#FF0000' \n";
}

void EpidemicAnalysis::mean_time_of_epidemic_over_lambda(std::vector<std::pair<double, double> > k_mean_time, int n) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_mean_time_by_lambda.gif";

    gp << "set title 'Tempo médio da epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Tempo médio' \n";
    gp << "set terminal gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'N=" + std::to_string(n) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_mean_time) << "title 'media' with linespoints lt rgb '#FF0000' \n";
}

void EpidemicAnalysis::std_time_of_epidemic_over_lambda(std::vector<std::pair<double, double> > k_std_time, int n) {
    Gnuplot gp;

    std::string filename = this->outputDir + "/epidemic_std_time_by_lambda.gif";

    gp << "set title 'Desvio padrão do tempo de epidemia com relação a k\n";
    gp << "set xlabel 'Número de passeios aleatórios' \n";
    gp << "set ylabel 'Desvio padrão' \n";
    gp << "set termina gif \n";
    std::cout << "filename: " << filename << "\n";
    gp << "set output '" + filename + "' \n";
    gp << "set key bottom right title 'N=" + std::to_string(n) + "' \n";
    //gp << "set logscale y \n";
    gp << "plot" << gp.file1d(k_std_time) << "title 'desvio padrao' with linespoints lt rgb '#FF0000' \n";
}