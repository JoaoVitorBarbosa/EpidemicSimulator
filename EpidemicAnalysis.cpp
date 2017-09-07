/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <set>
#include "EpidemicAnalysis.h"
#include <math.h>

EpidemicAnalysis::EpidemicAnalysis(std::string _outputDir) {
    outputDir = _outputDir;
}

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

void EpidemicAnalysis::analysisAll(std::string dir, int num_rw) {
    std::vector<std::pair<double, int> > all;

    for (int i = 0; i < num_rw; i++) {
        readTimestampStateChangeCSV(dir + "/RW_" + std::to_string(i) + "_Results.txt", all);
        all.push_back(std::make_pair(0, 0));
    }

    createStateStatistics(all, "System");
}

void EpidemicAnalysis::createStateStatistics(std::vector<std::pair<double, int> > timeStampStatesChange, std::string rw) {
    std::multiset<double> suscetibles_segments;
    std::multiset<double> contracteds_segments;
    std::multiset<double> infecteds_segments;

    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsInfected;
    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsContracted;
    boost::accumulators::accumulator_set<double, boost::accumulators::features < boost::accumulators::tag::mean, boost::accumulators::tag::median, boost::accumulators::tag::variance>> intervalsSusceptible;

    // states: 0 -> Susceptible, 1 -> Contracted, 2 -> infected
    // begins from 1 because first element is always time zero and initial state
    for (int i = 1; i < timeStampStatesChange.size(); i++) {
        // do not create negative interval for parse of all rw's
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

    std::string filename = this->outputDir + "/CCDFStateRW" + rw + ".png";
    
    // print all ccdf in one graphic
    Gnuplot gp;

    gp << "set title 'CCDF tempo gasto no estado \n";
    gp << "set xlabel 'Tempo gasto no estado' \n";
    gp << "set ylabel 'Fração' \n";
    //gp << "set yrange [" + std::to_string(1 / (double) size) + ":1]\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set logscale y \n";
    gp << "plot '-' title 'Susceptible' with lines lt rgb '#00FF00', "
            "'-' title 'Contracted' with lines lt rgb '#0000FF', "
            "'-' title 'Infected' with lines lt rgb '#FF0000' \n";

    gp.send1d(susceptiblesPts);
    gp.send1d(contractedsPts);
    gp.send1d(infectedsPts);

    std::ofstream arq;
    arq.open(this->outputDir + "/RW_" + rw + "_Results.txt");

    double mediaEsperada = 1/params.Rw.gama;
    double dpEsperado = mediaEsperada;
    double mediana = std::log(2)/params.Rw.gama;
    
    arq << "# Dados Infectado - Média: " << boost::accumulators::mean(intervalsInfected)
            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsInfected))
            << " Mediana: " << std::sqrt(boost::accumulators::median(intervalsInfected)) << '\n';

    arq << "# Dados Contraído - Média: " << boost::accumulators::mean(intervalsContracted)
            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsContracted))
            << " Mediana: " << std::sqrt(boost::accumulators::median(intervalsContracted)) << '\n';

    arq << "# Dados Suscetivel - Média: " << boost::accumulators::mean(intervalsSusceptible)
            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsSusceptible))
            << " Mediana: " << std::sqrt(boost::accumulators::median(intervalsSusceptible)) << '\n';

    arq.close();
}

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

void EpidemicAnalysis::printCCDF(std::multiset<double> elements, std::string type) {
    std::vector<std::pair<double, double> > xy_pts;
    int i = 0;
    int size = elements.size();
    for (std::multiset<double>::iterator it = elements.begin(); it != elements.end(); ++it) {
        double frac = (size - i) / (double) size;
        xy_pts.push_back(std::make_pair(*it, frac));
        i++;
    }

    Gnuplot gp;

    gp << "set title 'CCDF tempo gasto no estado " + type + "' \n";
    gp << "set xlabel 'Tempo gasto no estado' \n";
    gp << "set ylabel 'Fração' \n";
    gp << "set yrange [" + std::to_string(1 / (double) size) + ":1]\n";
    gp << "set term png \n";
    gp << "set output '../data/EpidemicAnalysis/" + type + ".png' \n";
    gp << "set logscale y \n";
    gp << "plot '-' with lines title '" + type + "'\n";

    gp.send1d(xy_pts);
}

void EpidemicAnalysis::randomWalkStateTimeSeries(std::string arquivo) {

    std::vector<std::pair<double, int> > infecteds;
    std::vector<std::pair<double, int> > contracteds;
    std::vector<std::pair<double, int> > susceptibles;

    std::ifstream arq;
    arq.open(arquivo.c_str());

    std::string linha = "";

    while (!arq.eof()) {

        std::string line;
        std::getline(arq, line);

        std::stringstream lineStream(line);
        std::string cell;

        if (line.find('#') == 0)
            continue;

        int i = 0;
        double time;
        int inf, sus, cont;
        while (std::getline(lineStream, cell, ',')) {
            if (i == 0)
                time = atof(cell.c_str());
            if (i == 1)
                inf = atof(cell.c_str());
            if (i == 2)
                cont = atof(cell.c_str());
            if (i == 3)
                sus = atof(cell.c_str());
            i++;
        }

        infecteds.push_back(std::make_pair(time, inf));
        contracteds.push_back(std::make_pair(time, cont));
        susceptibles.push_back(std::make_pair(time, sus));
    }

    arq.close();

    std::string filename = this->outputDir + "/rwStateTimeSeries.png";
    
    Gnuplot gp;

    gp << "set title '' \n";
    gp << "set xlabel 'Time' \n";
    gp << "set ylabel 'Number of RW' \n";
    //gp << "set xrange [0:" + std::to_string(max + 10) + "]\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    gp << "set terminal png medium size 1400,900 \n";
    gp << "set yrange [0:" + std::to_string(params.Rw.k) + "]\n";
    gp << "plot '-' title 'Susceptible' with lines lt rgb '#00FF00', '-' title 'Contracted' with lines lt rgb '#0000FF', '-' title 'Infected' with lines lt rgb '#FF0000' \n";

    gp.send1d(susceptibles);
    gp.send1d(contracteds);
    gp.send1d(infecteds);

}

void VertexTimeInfecteds(int vertex) {
    std::string v = std::to_string(vertex);
    std::vector<std::pair<double, int> > pairs;
    //readCSV("../data/output/" + v +"_vertexResults.txt", pairs);
}

void readTimeNumberInfect(std::string filename, std::set<double, int>) {
    std::ifstream arq;
    arq.open(filename.c_str());

    std::string linha = "";

    //std::set 

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

void EpidemicAnalysis::analysisStateTime(std::string filepath, std::string title) {
    std::vector<std::pair<double, int> > pairs;
    readTimestampStateChangeCSV(filepath, pairs);
    createStateStatistics(pairs, title);
}

void EpidemicAnalysis::infectedsGraphic(std::string filepath, std::string title) {
    std::vector<std::pair<int,double>> timeStampStatesChange;
    
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
    
    std::string filename = outputDir + "/infecteds" + title + ".png";
    //../data/EpidemicAnalysis/infectInterval.png
    
    Gnuplot gp;

    gp << "set title 'Fração do tempo do sistema com k infectados' \n";
    gp << "set xlabel 'Numero de infectados' \n";
    gp << "set ylabel 'Fração' \n";
    //gp << "set yrange [" + std::to_string(1 / (double) size) + ":1]\n";
    gp << "set term png \n";
    gp << "set output '" + filename + "' \n";
    //gp << "set logscale y \n";
    gp << "plot '-' with lines title 'infectados'\n";

    gp.send1d(timeStampStatesChange);   
    
}