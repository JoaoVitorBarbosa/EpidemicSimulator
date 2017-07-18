/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "Simulator.h"
#include <unistd.h>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <ctime>
#include <map>

#define GetCurrentDir getcwd

std::string GetCurrentWorkingDir(void) {
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}

void Simulator::changeNumberInfected(int num_Infec) {
    double interval = time - timeLastNumberInfect;
    double actInterval = infectedInterval[num_Infecteds];
    infectedInterval[num_Infecteds] = interval + actInterval;

    timeLastNumberInfect = time;
    num_Infecteds = num_Infec;
}

std::string Simulator::eventToString(EventType evt) {
    switch (evt) {
        case EventType::Infect:
            return "Infect";
            break;
        case EventType::Recover:
            return "Recover";
            break;
        case EventType::Walk:
            return "Walk";
            break;
    }
}

void Simulator::fillVertices(Params params) {
    for (int i = 0; i < graph.num_vertices; i++) {
        Vertex * vertex = new Vertex(0, params.p, i);
        vertices.push_back(vertex);
    }
}

void Simulator::setupRandomWalks(RwParam rwParams) {
    std::map<int, bool> rwToInf;
    for (int j = 0; j < rwParams.ki; j++) {
        int rw = rg.uniform(0, rwParams.ki);
        while (!(rwToInf[rw] == NULL))
            rw = rg.uniform(0, rwParams.ki);
        rwToInf[rw] = true;
    }

    for (int i = 0; i < k; i++) {
        std::vector<double> rwParameters = rwParams.rwParamVector[i];
        int vp = rwParameters.empty() || rwParameters[1] == -1 ? rg.uniform(0, graph.num_vertices) : rwParameters[1]; 
        //RandomWalk * rw = new RandomWalk(vp, 1, 0.05, .5, s, i);
        State s = rwToInf[i] ? State::Infected : State::Susceptible;
        RandomWalk * rw = !rwParameters.empty() ?
                new RandomWalk(vp, rwParameters[2], rwParameters[3], rwParameters[4], (State) rwParameters[5], (int) rwParameters[0], outputDir)
                : new RandomWalk(vp, rwParams.lambda, rwParams.gama, rwParams.tau, s, i, outputDir);
        randomWalks.push_back(rw);

        rw->setRwPosition(vertices[vp]->setRandomWalk(rw));
        if (rw->getState() == State::Infected) {
            double t_recover = time + rg.exponential(rw->getGama());
            Event evtRec(t_recover, i, EventType::Recover);
            events.push(evtRec);

            rw->insertTimeInfected(t_recover - time);

            vertices[vp]->increaseRwInfecteds(time);
            changeNumberInfected(num_Infecteds + 1);

            rw->setTimeStateChange(0.0, 2);
        } else
            rw->setTimeStateChange(0.0, 0);

        //walk event
        double timeToWalk = rg.exponential(randomWalks[i]->getLambda());
        rw->insertTimeWalking(timeToWalk);

        double walkTime = time + timeToWalk;
        Event evt(walkTime, i, EventType::Walk);
        events.push(evt);

        // write first event (creation)
        History history(rw->getCode(), rw->getVertex(), time, rw->stateToString(), eventToString(evt.type), "-");
        rw->writeEvent(history);
    }

}

void Simulator::fillEvents() {
    for (int i = 0; i < randomWalks.size(); i++) {
        double walkTime = rg.exponential(randomWalks[i]->getLambda());
        Event evt(walkTime, i, EventType::Walk);

        events.push(evt);
    }
}

Simulator::Simulator(Params params, std::string jsonStr, ManipulaGrafoV _graph) {

    k = params.Rw.k;
    timeLimit = params.Time;
    rounds = params.Rounds;
    // implementar Graph -> Inserir direto Classe ManipulaGrafo com os vertices e arestas
    this->graph = _graph;
    time = 0;
    num_Inf_Events = 0;
    num_Infecteds = 0;
    timeLastNumberInfect = 0;

    // create new directory for each execution
    time_t now = std::time(0);
    tm *ltm = localtime(&now);
    std::string fName = std::to_string(ltm->tm_mday) + "-" + std::to_string(1 + ltm->tm_mon) + "-" + std::to_string(1900 + ltm->tm_year) + " " + std::to_string(ltm->tm_hour) + ":" + std::to_string(ltm->tm_min) + ":" + std::to_string(1 + ltm->tm_sec);

    outputDir = params.OutputDir + "/" + fName;
    boost::filesystem::path dir(outputDir.c_str());
    if (boost::filesystem::create_directory(dir)) {
        std::cout << "Directory created successfully" << "\n";
    }

    fileNameInfectInterval = outputDir + "/infectedsInterval.txt";
    
    // create file to store results
    fileNameNumberRandomWalkStates = outputDir + "/system_results.txt";
    std::ofstream arq;
    arq.open(fileNameNumberRandomWalkStates);
    arq << "# time,infected,contracted,susceptible" << "\n";
    arq.close();

    // writting params to file
    std::string paramsFile = outputDir + "/Params_Execution.txt";
    arq.open(paramsFile);
    arq << jsonStr << "\n";
    arq.close();

    // initialize parameters
    initialize(params);
}

Simulator::Simulator(int _k, int _timeLimt, int _rounds, std::string rede) {
    k = _k;
    timeLimit = _timeLimt;
    rounds = _rounds;
    graph.lerArquivo(rede, false);
    time = 0;
    timeLastNumberInfect = 0;
    num_Inf_Events = 0;
    num_Infecteds = 0;

    fileNameNumberRandomWalkStates = "../data/output/system_results.txt";
    std::ofstream arq;
    arq.open(fileNameNumberRandomWalkStates);
    arq << "# time,infected,contracted,susceptible" << "\n";
    arq.close();
}

Simulator::~Simulator() {
    for (std::vector<Vertex*>::iterator it = vertices.begin(); it != vertices.end(); it++)
        delete *it;

    for (std::vector<RandomWalk*>::iterator it = randomWalks.begin(); it != randomWalks.end(); it++)
        delete *it;
}

Event Simulator::getEventTop() {
    Event evt = events.top();
    events.pop();
    return evt;
}

void Simulator::initialize(Params params) {
    // fill vertices
    fillVertices(params);
    // fill rws and events
    setupRandomWalks(params.Rw);
    //fill Events
    //fillEvents();
}

void Simulator::infect(Vertex * vertex, Event evt) {
    // all randomwalk are infecteds
    if (vertex->getRwInfecteds() == vertex->getRandomWalkList().size())
        return;

    std::list<RandomWalk*>::iterator it = vertex->getRandomWalkList().begin();
    for (int i = 0; i < vertex->getRandomWalkList().size(); i++) {
        if ((*it)->getState() == State::Susceptible) {
            bool bern = rg.bernoulli(vertex->getP());
            if (bern) {
                double t_tau = time + rg.exponential((*it)->getTau());
                Event evt_new(t_tau, (*it)->getCode(), EventType::Infect);
                events.push(evt_new);

                (*it)->setState(State::Contracted);
                (*it)->setTimeStateChange(time, State::Contracted);
                num_Inf_Events++;

                (*it)->insertTimeContracted(t_tau - time);

                History history((*it)->getCode(), (*it)->getVertex(), time, (*it)->stateToString(), "-", "Infected", evt.randomwalk);
                (*it)->writeEvent(history);
            }
        }
        it++;
    }

}

void Simulator::beInfected(Vertex * v, Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    bool inf = 1 - std::pow((1 - v->getP()), v->getRwInfecteds());
    if (inf) {
        double t_tau = time + rg.exponential(rw->getTau());
        Event evt_new(t_tau, evt.randomwalk, EventType::Infect);
        events.push(evt_new);

        rw->setState(State::Contracted);
        rw->setTimeStateChange(time, State::Contracted);
        num_Inf_Events++;

        rw->insertTimeContracted(t_tau - time);

        History history(rw->getCode(), rw->getVertex(), time, rw->stateToString(), "-", "Infected");
        rw->writeEvent(history);
    }
}

void Simulator::processWalk(Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);

    double timeToWalk = rg.exponential(rw->getLambda());
    rw->insertTimeWalking(timeToWalk);

    double t_walk = timeToWalk + time;

    // create new event and insert in events pq
    Event evt_new(t_walk, evt.randomwalk, EventType::Walk);
    events.push(evt_new);

    // get neighbor to go
    int vIndex = rw->getVertex();
    Vertex * srcVertex = vertices.at(vIndex);
    std::vector<unsigned int> neighbors = graph.get_vizinhos(vIndex + 1);
    int degree = graph.getGrau(vIndex + 1);
    int indexN = rg.uniform(0, degree); //uniform will not return b parameter
    int vertexNeighbor = neighbors.at(indexN);

    // delete rw in old vertex
    srcVertex->eraseRandomWalk(rw->getRwPosition());

    // vertex needs to be a pointer, if not gets a copy and don't change the real vertex in vector
    Vertex * dstVertex = vertices.at(vertexNeighbor - 1);
    randomWalks.at(evt.randomwalk)->setRwPosition(dstVertex->setRandomWalk(rw));
    randomWalks.at(evt.randomwalk)->setVertex(vertexNeighbor - 1);

    //std::cout << "Vertex to go: " << vertexNeighbor -1 << std::endl;

    History history(rw->getCode(), rw->getVertex(), time, rw->stateToString(), eventToString(evt.type), "-");
    rw->writeEvent(history);


    if (rw->getState() == State::Infected) {
        srcVertex->decreaseRwInfecteds(time);
        dstVertex->increaseRwInfecteds(time);
        infect(dstVertex, evt);
    } else if (rw->getState() == State::Susceptible)
        beInfected(dstVertex, evt);
}

void Simulator::processInfect(Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    rw->setState(State::Infected);
    rw->setTimeStateChange(time, State::Infected);

    // generate time to recover
    double t_recover = time + rg.exponential(rw->getGama());
    Event evt_new(t_recover, evt.randomwalk, EventType::Recover);
    events.push(evt_new);

    History history(rw->getCode(), rw->getVertex(), time, rw->stateToString(), eventToString(evt.type), "-");
    rw->writeEvent(history);
    rw->insertTimeInfected(t_recover - time);

    // increase number of vertices infected
    Vertex * vertex = vertices.at(rw->getVertex());
    vertex->increaseRwInfecteds(time);

    changeNumberInfected(num_Infecteds + 1);
    num_Inf_Events--;

    infect(vertex, evt);
}

void Simulator::processHeal(Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    rw->setState(State::Susceptible);
    rw->setTimeStateChange(time, State::Susceptible);

    History history(rw->getCode(), rw->getVertex(), time, rw->stateToString(), eventToString(evt.type), "-");
    rw->writeEvent(history);

    changeNumberInfected(num_Infecteds - 1);

    Vertex * vertex = vertices.at(rw->getVertex());
    vertex->decreaseRwInfecteds(time);

    beInfected(vertex, evt);
}

void Simulator::writeNumberRwStatePerTime() {
    std::ofstream arq;
    arq.open(fileNameNumberRandomWalkStates, std::ofstream::out | std::ofstream::app);

    arq << std::fixed << time << "," << num_Infecteds << "," << num_Inf_Events << "," << k - (num_Infecteds + num_Inf_Events) << std::endl;

    arq.close();
}

void Simulator::writeInfectInterval(){
    std::ofstream arq;
    arq.open(fileNameInfectInterval, std::ofstream::out | std::ofstream::app);

    for(int i = 0; i < k; i++)
        arq << i << "," << infectedInterval[i]/time << std::endl;

    arq.close();
}

void Simulator::process() {
    std::cout << "Begin epidemic" << std::endl;
    double old_percentual = 0.0;
    while (!events.empty()) {

        if (time > timeLimit || (num_Inf_Events == 0 && num_Infecteds == 0))
            break;

        Event evt = getEventTop();
        time = evt.time;

        switch (evt.type) {
            case Walk:
                processWalk(evt);
                break;
            case Infect:
                processInfect(evt);
                break;
            case Recover:
                processHeal(evt);
                break;

            default:
                break;
        }

        if (evt.type != EventType::Walk)
            writeNumberRwStatePerTime();

        double percentual = time / timeLimit * 100;

        if ((int) percentual != (int) old_percentual)
            std::cout << std::fixed << percentual << "%" << std::endl;

        old_percentual = percentual;
    }

    std::cout << "End epidemic" << std::endl;

    for (int k = 0; k < randomWalks.size(); k++) {
        randomWalks.at(k)->WriteResults();
    }
    
    writeInfectInterval();
    
    /*while(true){
    
        std::cout << "Digite o numero do passeio." << std::endl;
        std::cin >> stop;
        
        if(stop == -1)
            break;
        
        if(stop > k)
            continue;
        
        randomWalks.at(stop)->drawStates();
    
        std::cin.get();
    }
    
    std::cout << std::endl;*/


}