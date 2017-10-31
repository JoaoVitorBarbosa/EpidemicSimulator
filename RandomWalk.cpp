/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "RandomWalk.h"
#include "Simulator.h"

RandomWalk::RandomWalk() {
}

RandomWalk::RandomWalk(int _vertex, double _lambda, double _gama, double _tau, State _state, int _code, std::string _outputDir) {
    vertex = _vertex;
    lambda = _lambda;
    gama = _gama;
    tau = _tau;
    state = _state;
    code = _code;
    outputDir = _outputDir;

    std::string prefix = "RW" + std::to_string(code);
    fileName = outputDir + "/" + prefix + "_change  .txt";
    fileNameParmResult = outputDir + "/" + prefix + "_Results.txt";
    fileNameWalkingTimes = outputDir + "/" + prefix + "_walkingTimes.txt";

    std::ofstream arq;

    arq.open(fileName);
    arq << "# Random Walk: " << code << std::endl;
    arq << "# Parameters: Lambda: " << lambda << ", Gama: " << gama << ", Tau: " << tau << '\n';
    arq << "# time,event,vertex,state,effect,infected_by" << '\n';

    arq.close();

    arq.open(fileNameParmResult);
    arq << "# Parameters: Lambda: " << lambda << ", Gama: " << gama << ", Tau: " << tau << '\n';
    arq << "# time,state_changed" << '\n';
    arq.close();
    
    arq.open(fileNameWalkingTimes);
    arq << "# Walking Times" << std::endl;
    arq.close();


}

void RandomWalk::setVertex(int v) {
    vertex = v;
}

int RandomWalk::getVertex() {
    return vertex;
}

void RandomWalk::setCode(int i) {
    code = i;
}

int RandomWalk::getCode() {
    return code;
}

void RandomWalk::setState(State _state) {
    state = _state;
}

State RandomWalk::getState() {
    return state;
}

void RandomWalk::setLambda(double _lambda) {
    lambda = _lambda;
}

double RandomWalk::getLambda() {
    return lambda;
}

void RandomWalk::setGama(double _gama) {
    gama = _gama;
}

double RandomWalk::getGama() {
    return gama;
}

void RandomWalk::setTau(double _tau) {
    tau = _tau;
}

double RandomWalk::getTau() {
    return tau;
}

void RandomWalk::setRwPosition(std::list<RandomWalk*>::iterator it) {
    rwIt = it;
}

std::list<RandomWalk*>::iterator RandomWalk::getRwPosition() {
    return rwIt;
}

bool RandomWalk::isInfected() {
    return state == Contracted;
}

std::string RandomWalk::stateToString() {
    switch (state) {
        case State::Contracted:
            return "Contracted";
        case State::Infected:
            return "Infected";
        case State::Susceptible:
            return "Susceptible";
    }
}

void RandomWalk::writeEvent(History history) {
    //std::cout << "RW: " << code << " ; Vertex: " << vertex << " ; Time: " << time << " ; Event: " << event.toString() << " ; State: " << stateToString() << std::endl;
    //writeFile(event, time);
    writeFile(history);
}

void RandomWalk::writeFile(History history) {
//    std::ofstream arq;
//    arq.open(fileName, std::ofstream::out | std::ofstream::app);
//    arq << std::fixed << history.time << "," << history.event << "," << history.vertex << "," << history.state << "," << history.effect << "," << history.getRWInfected() << '\n';
//    arq.close();
//    
//    if(history.event == "Walk")
//    {
//        arq.open(fileNameWalkingTimes, std::ofstream::out | std::ofstream::app);
//        arq << std::fixed << history.time << std::endl;
//        arq.close();
//    }
}

void RandomWalk::insertTimeWalking(double _t) {
    intervalsWalking(_t);
}

void RandomWalk::insertTimeInfected(double _t) {
    intervalsInfected(_t);
}

void RandomWalk::insertTimeContracted(double _t) {
    intervalsContracted(_t);
}

void RandomWalk::WriteResults() {
//    std::ofstream arq;
//
//    arq.open(fileNameParmResult, std::ofstream::out | std::ofstream::app);
//
//    arq << std::fixed << "# Dados Infectado - Média: " << boost::accumulators::mean(intervalsInfected)
//            << " Mediana: " << boost::accumulators::median(intervalsInfected) 
//            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsInfected)) << '\n';
//    
//    arq << std::fixed << "# Dados Contraído - Média: " << boost::accumulators::mean(intervalsContracted) 
//            << " Mediana: " << boost::accumulators::median(intervalsContracted) 
//            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsContracted)) << '\n';
//    
//    arq << std::fixed << "# Dados Caminhar - Média: " << boost::accumulators::mean(intervalsWalking) 
//            << " Mediana: " << boost::accumulators::median(intervalsWalking) 
//            << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsWalking)) << '\n';
//
//    arq.close();
}

void RandomWalk::setTimeStateChange(double _t, int _s) {
//    timeStampStatesChange.push_back(std::make_pair(_t, _s));
//
//    std::ofstream arq;
//    arq.open(fileNameParmResult, std::ofstream::out | std::ofstream::app);
//
//    arq << std::fixed << _t << "," << _s << '\n';
//
//    arq.close();
}

void RandomWalk::drawStates() {
    Gnuplot gp;

    std::vector<std::vector<std::pair<double, int> > > suscetibles_segments;
    std::vector<std::vector<std::pair<double, int> > > contracteds_segments;
    std::vector<std::vector<std::pair<double, int> > > infecteds_segments;

    for (int i = 1; i < timeStampStatesChange.size(); i++) {
        // time infected
        if (timeStampStatesChange[i].second == State::Susceptible) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(timeStampStatesChange[i - 1].first, 3));
            segment.push_back(std::make_pair(timeStampStatesChange[i].first, 3));
            infecteds_segments.push_back(segment);
        }

        // time suscetible
        if (timeStampStatesChange[i].second == State::Contracted) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(timeStampStatesChange[i - 1].first, 1));
            segment.push_back(std::make_pair(timeStampStatesChange[i].first, 1));
            suscetibles_segments.push_back(segment);
        }

        // time contracted
        if (timeStampStatesChange[i].second == State::Infected) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(timeStampStatesChange[i - 1].first, 2));
            segment.push_back(std::make_pair(timeStampStatesChange[i].first, 2));
            contracteds_segments.push_back(segment);
        }
    }

    int max = timeStampStatesChange[timeStampStatesChange.size() - 1].first;

    gp << "set title 'Tempo do RW " + std::to_string(code) + " em cada Estado' \n";
    gp << "set xlabel 'Time' \n";
    gp << "set ylabel 'State \n'";
    gp << "set xrange [0:" + std::to_string(max + 10) + "]\n";
    gp << "set yrange [0:4]\n";
    gp << "plot '-' title 'Susceptible' with linespoints lt rgb '#00FF00', '-' title 'Contracted' with linespoints lt rgb '#0000FF', '-' title 'Infected' with linespoints lt rgb '#FF0000' \n";
    //gp << "plot '-' with linespoints lt rgb '#FF0000', '-' with linespoints lt rgb '#00FF00'\n";

    gp.send2d(suscetibles_segments);
    gp.send2d(contracteds_segments);
    gp.send2d(infecteds_segments);

    //	std::cin.get();


}