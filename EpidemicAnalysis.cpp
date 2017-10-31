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
    do_validation = false;
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

    std::string prefix("RW" + rw);
    std::string filename = this->outputDir + "/" + prefix + "_CCDFState.png";

    // print all ccdf in one graphic
    Gnuplot gp;

    if (do_validation) {
        gp << "set title 'CCDF tempo gasto no estado infectado' \n";
        gp << "set xlabel 'Tempo' \n";
        gp << "set ylabel 'Fração' \n";
        gp << "set term png \n";
        gp << "set output '" + this->outputDir + "/" + prefix + "_CCDFInfected" + "' \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(infectedsPts) << "title 'Infectado' with lines lt rgb '#00FF00', exp(-" + std::to_string(gama) + "*x) title 'Statistics' with lines lt rgb '#FF0000' \n";

        gp << "set title 'CCDF tempo gasto no estado contraido' \n";
        gp << "set xlabel 'Tempo' \n";
        gp << "set ylabel 'Fração' \n";
        gp << "set term png \n";
        gp << "set output '" + this->outputDir + "/" + prefix + "_CCDFContracted" + "' \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(contractedsPts) << "title 'Contraido' with lines lt rgb '#00FF00', exp(-" + std::to_string(tau) + "*x) title 'Statistics' with lines lt rgb '#FF0000' \n";
    } else {
        gp << "set title 'CCDF tempo gasto no estado \n";
        gp << "set xlabel 'Tempo gasto no estado' \n";
        gp << "set ylabel 'Fração' \n";
        gp << "set term png \n";
        gp << "set output '" + filename + "' \n";
        gp << "set logscale y \n";
        gp << "plot" << gp.file1d(susceptiblesPts) << "title 'Susceptible' with lines lt rgb '#00FF00', "
                << gp.file1d(contractedsPts) << "title 'Contracted' with lines lt rgb '#0000FF', "
                << gp.file1d(infectedsPts) << "title 'Infected' with lines lt rgb '#FF0000' \n";
    }

    std::ofstream arq;
    arq.open(this->outputDir + "/" + prefix + "_Results.txt", std::ofstream::out | std::ofstream::app);

    double infected_mean = boost::accumulators::mean(intervalsInfected);

    arq << "Estado Infectado \n"
            << "    Média teórica:          " << 1 / gama << "    Média estatística: " << infected_mean << "\n"
            << "    Desvio Padrão teórico:  " << 1 / gama << "    Desvio padrão estatístico: " << std::sqrt(boost::accumulators::variance(intervalsInfected)) << "\n"
            << "    Mediana teórica:        " << std::log(2) / gama << "    Mediana estatística: " << get_median(infecteds_segments).median << "\n"
            << "    Gama teórico:           " << gama << "    Gama estatístico: " << 1 / infected_mean << "\n"
            << '\n';

    arq << "Estado Contraído \n"
            << "    Média teórica:          " << 1 / tau << "    Média estatística: " << boost::accumulators::mean(intervalsContracted) << "\n"
            << "    Desvio Padrão teórico:  " << 1 / tau << "    Desvio padrão estatístico: " << std::sqrt(boost::accumulators::variance(intervalsContracted)) << "\n"
            << "    Mediana teórica:        " << std::log(2) / tau << "    Mediana estatística: " << get_median(contracteds_segments).median << "\n"
            << "    Tau teórico:           " << tau << "    Tau estatístico: " << 1 / boost::accumulators::mean(intervalsContracted) << "\n"
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

    std::string filename = this->outputDir + "/" + "RW" + randomWalk + "_CCDF_walk.png";

    Gnuplot gp;

    gp << "set title 'CCDF do tempo de caminhada\n";
    gp << "set xlabel 'Tempo até próximo passo' \n";
    gp << "set ylabel 'Fração' \n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set logscale y \n";
    if (do_validation) {
        gp << "plot" << gp.file1d(walkingPts) << "title 'empiric' with lines lt rgb '#00FF00', exp(-" + std::to_string(lambda) + "*x) title 'statistics' with lines lt rgb '#FF0000' \n";

        // writing in file
        std::string prefix("RW" + randomWalk);
        std::string filename = this->outputDir + "/" + prefix + "_CCDFState.png";

        std::ofstream arq;
        arq.open(this->outputDir + "/" + prefix + "_Results.txt", std::ofstream::out | std::ofstream::app);

        double emp_mean = boost::accumulators::mean(intervalsWalking);
        Stats stats = get_median(walking_segments);

        arq << "Dados do evento caminhar \n"
                << "    Média teórica:          " << 1 / lambda << "    Média estatística: " << emp_mean << "   Outra média: " << stats.mean << "\n"
                << "    Desvio Padrão teórico:  " << 1 / lambda << "    Desvio padrão estatístico: " << std::sqrt(boost::accumulators::variance(intervalsWalking)) << "\n"
                << "    Mediana teórica:        " << std::log(2) / lambda << "    Mediana estatística: " << stats.median << "\n"
                << "    Lambda teórico:           " << lambda << "    Gama estatístico: " << 1 / emp_mean << "\n"
                << '\n';
    } else
        gp << "plot" << gp.file1d(walkingPts) << "title 'empiric' with lines lt rgb '#00FF00 \n'";

}

void EpidemicAnalysis::randomWalkStateTimeSeries(std::string arquivo) {

    std::vector<std::pair<double, int> > infecteds;
    std::vector<std::pair<double, int> > contracteds;
    std::vector<std::pair<double, int> > susceptibles;

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

        infecteds.push_back(std::make_pair(time, inf));
        contracteds.push_back(std::make_pair(time, cont));
        susceptibles.push_back(std::make_pair(time, sus));
    }

    arq.close();

    max = max + 1;

    std::string filename = this->outputDir + "/rwStateTimeSeries.png";

    Gnuplot gp;
    gp << "set title '' \n";
    gp << "set xlabel 'Time' \n";
    gp << "set ylabel 'Number of RW' \n";
    //gp << "set xrange [0:" + std::to_string(max + 10) + "]\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set terminal png medium size 1400,900 \n";
    gp << "set yrange [0:" + std::to_string(max) + "]\n";
    gp << "plot '-' title 'Susceptible' with lines lt rgb '#00FF00', '-' title 'Contracted' with lines lt rgb '#0000FF', '-' title 'Infected' with lines lt rgb '#FF0000' \n";

    gp.send1d(susceptibles);
    gp.send1d(contracteds);
    gp.send1d(infecteds);

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

    gp << "set title 'Fração do tempo do sistema com k infectados' \n";
    gp << "set xlabel 'Numero de infectados' \n";
    gp << "set ylabel 'Fração' \n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    //gp << "set logscale y \n";
    gp << "plot '-' title 'Susceptible' with lines lt rgb '#00FF00', '-' title 'Contracted' with lines lt rgb '#0000FF', '-' title 'Infected' with lines lt rgb '#FF0000' \n";

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