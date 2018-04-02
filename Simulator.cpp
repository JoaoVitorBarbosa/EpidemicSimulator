/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "Simulator.h"
#include "Logger.h"

#define GetCurrentDir getcwd

Simulator::Simulator(Params params, std::string jsonStr, ManipulaGrafoV _graph) {
    do_analysis = true;
    k = params.Rw.k;
    limit_time_epidemic = params.Time;
    rounds = params.Runs;
    continue_simulation_after_epidemic = params.cont_simulation;
    // implementar Graph -> Inserir direto Classe ManipulaGrafo com os vertices e arestas
    this->graph = _graph;

    time = 0;
    num_Inf_Events = 0;
    num_Infected = 0;
    num_Contracted = 0;
    time_of_last_number_of_infected = 0;
    time_of_last_number_of_contracted = 0;
    time_of_last_number_of_susceptible = 0;

    outputDir = params.OutputDir; // + "/" + fName;
    boost::filesystem::path dir(outputDir.c_str());
    if (boost::filesystem::create_directory(dir)) {
        Logger::Info("Directory created successfully\n");
    }

    fileNameInfectInterval = outputDir + "/infectedsInterval.txt";
    fileNameContractedInterval = outputDir + "/contractedInterval.txt";
    fileNameSusceptibleInterval = outputDir + "/susceptibleInterval.txt";
    file_name_system_results = outputDir + "/epidemic_results.txt";

    // create file to store results
    fileNameNumberRandomWalkStates = outputDir + "/system_results.txt";
    std::ofstream arq;
    arq.open(fileNameNumberRandomWalkStates);
    arq << "# time,infected,contracted,susceptible" << "\n";
    arq.close();

    // writting params to file
    std::string paramsFile = outputDir + "/Params_Execution.json";
    arq.open(paramsFile);
    arq << jsonStr << "\n";
    arq.close();

    // file to write infected density
    file_infected_density = outputDir + "/infected_density.csv";
    arq.open(file_infected_density);
    arq << "#time,infected_density" << "\n";
    arq.close();

    // initialize parameters
    initialize(params);
}

Simulator::Simulator(int _k, int _timeLimt, int _rounds, std::string rede) {
    k = _k;
    limit_time_epidemic = _timeLimt;
    rounds = _rounds;
    graph.lerArquivo(rede, false);
    time = 0;
    time_of_last_number_of_infected = 0;
    time_of_last_number_of_contracted = 0;
    time_of_last_number_of_susceptible = 0;
    num_Inf_Events = 0;
    num_Infected = 0;

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

std::string GetCurrentWorkingDir(void) {
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}

void Simulator::change_infected_number(int num_Infec) {
    double interval = time - time_of_last_number_of_infected;
    double cInterval = infected_intervals[num_Infected];
    infected_intervals[num_Infected] = interval + cInterval;

    time_of_last_number_of_infected = time;

    // number of infected decreases and susceptible increase
    if (num_Infec < num_Infected) {
        interval = time - time_of_last_number_of_susceptible;
        int num_Susceptible = k - num_Infected - num_Contracted;
        cInterval = susceptible_intervals[num_Susceptible];
        susceptible_intervals[num_Susceptible] = interval + cInterval;

        time_of_last_number_of_susceptible = time;
    }

    num_Infected = num_Infec;
}

void Simulator::changeNumberContracted(int num_Cont) {
    double interval = time - time_of_last_number_of_contracted;
    double cInterval = contracted_intervals[num_Contracted];
    contracted_intervals[num_Contracted] = interval + cInterval;

    time_of_last_number_of_contracted = time;

    // number of infected decreases and susceptible increase
    if (num_Cont > num_Contracted) {
        interval = time - time_of_last_number_of_susceptible;
        int num_Susceptible = k - num_Infected - num_Contracted;
        cInterval = susceptible_intervals[num_Susceptible];
        susceptible_intervals[num_Susceptible] = interval + cInterval;

        time_of_last_number_of_susceptible = time;
    }

    num_Contracted = num_Cont;
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

void Simulator::fill_vertices(Params params) {
    // create directory to vertices 
    // remove: don't create vertices directory to save space and be faster
    std::string verticesOutputDir = this->outputDir + "/Vertices";
    if (do_analysis)
        Utils::createDirectory(verticesOutputDir);

    for (int i = 0; i < graph.num_vertices; i++) {
        double p = params.Vertex.vertexParamVector[i] != NULL ? params.Vertex.vertexParamVector[i] : params.Vertex.p;
        Vertex * vertex = new Vertex(0, p, i, k, verticesOutputDir);
        vertices.push_back(vertex);
    }
}

void Simulator::setupRandomWalks(RwParam rwParams) {
    // create directory to random walks
    // remove after: disable writing for random walk to save space and be faster
    std::string randomWalkOutputDir = this->outputDir + "/RandomWalks";
    if (do_analysis)
        Utils::createDirectory(randomWalkOutputDir);

    // select randomly random walks that must be infected
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
        State s = rwToInf[i] ? State::Infected : State::Susceptible;
        RandomWalk * rw = !rwParameters.empty() ?
                new RandomWalk(vp, rwParameters[2], rwParameters[3], rwParameters[4], (State) rwParameters[5], (int) rwParameters[0], randomWalkOutputDir)
                : new RandomWalk(vp, rwParams.lambda, rwParams.gama, rwParams.tau, s, i, randomWalkOutputDir);
        randomWalks.push_back(rw);

        rw->set_rw_position(vertices[vp]->setRandomWalk(rw));

        if (rw->get_state() == State::Infected) {
            double t_recover = time + rg.exponential(rw->get_gama());
            Event evtRec(t_recover, i, EventType::Recover);
            events.push(evtRec);

            rw->insertTimeInfected(t_recover - time);

            vertices[vp]->increaseRwInfecteds(time);
            change_infected_number(num_Infected + 1);

            rw->setTimeStateChange(0.0, 2);
        } else
            rw->setTimeStateChange(0.0, 0);

        if (rw->get_state() == State::Contracted)
            changeNumberContracted(num_Contracted + 1);

        //walk event
        double timeToWalk = rg.exponential(randomWalks[i]->get_lambda());
        rw->insertTimeWalking(timeToWalk);

        double walkTime = time + timeToWalk;
        Event evt(walkTime, i, EventType::Walk);
        events.push(evt);

        // write first event (creation)
        History history(rw->get_code(), rw->get_vertex(), time, rw->state_to_string(), eventToString(evt.type), "-");
        rw->writeEvent(history);
    }

}

void Simulator::fillEvents() {
    for (int i = 0; i < randomWalks.size(); i++) {
        double walkTime = rg.exponential(randomWalks[i]->get_lambda());
        Event evt(walkTime, i, EventType::Walk);

        events.push(evt);
    }
}

Event Simulator::get_top_event() {
    Event evt = events.top();
    events.pop();
    return evt;
}

void Simulator::initialize(Params params) {
    // fill vertices
    fill_vertices(params);
    // fill rws and events
    setupRandomWalks(params.Rw);
    //fill Events
    //fillEvents();
}

void Simulator::infect(Vertex * vertex, Event evt) {
    // if all randomwalks are infected, leave

    //Logger::Trace("Infected: " + std::to_string(vertex->getRwInfecteds()) + " Total: " + std::to_string(vertex->getRandomWalkList().size()));

    if (vertex->getRwInfecteds() == vertex->getRandomWalkList().size())
        return;

    //Logger::Trace("For of infect: " + std::to_string(vertex->getRandomWalkList().size()));

    std::list<RandomWalk*>::iterator it = vertex->getRandomWalkList().begin();
    for (int i = 0; i < vertex->getRandomWalkList().size(); i++) {
        if ((*it)->get_state() == State::Susceptible) {
            // increase encounters number
            vertex->increase_encounters_by(1);

            bool bern = rg.bernoulli(vertex->getP());
            if (bern) {

                // encounter generated infection
                vertex->increase_encounters_with_transmition_by(1);

                double t_tau = time + rg.exponential((*it)->get_tau());
                Event evt_new(t_tau, (*it)->get_code(), EventType::Infect);
                events.push(evt_new);

                // write snapshot of states
                writeNumberRwStatePerTime("contract"); //evt.toString());

                (*it)->set_state(State::Contracted);
                (*it)->setTimeStateChange(time, State::Contracted);
                num_Inf_Events++;
                infectionTimes++; // increase number of infections

                (*it)->insertTimeContracted(t_tau - time);

                History history((*it)->get_code(), (*it)->get_vertex(), time, (*it)->state_to_string(), "-", "Infected", evt.randomwalk);
                (*it)->writeEvent(history);

                // increase number of contracted
                changeNumberContracted(num_Contracted + 1);
            }
        }
        it++;
    }

}

void Simulator::beInfected(Vertex * v, Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    double p = v->getP();
    int infected = v->getRwInfecteds();
    double new_p = 1 - std::pow((1 - p), infected);
    bool inf = rg.bernoulli(new_p);

    if (inf) {
        v->sum_success_encounters(infected);

        double t_tau = time + rg.exponential(rw->get_tau());
        Event evt_new(t_tau, evt.randomwalk, EventType::Infect);
        events.push(evt_new);

        // write snapshot of states
        writeNumberRwStatePerTime("contract"); //evt.toString());

        rw->set_state(State::Contracted);
        rw->setTimeStateChange(time, State::Contracted);
        num_Inf_Events++;
        infectionTimes++;

        rw->insertTimeContracted(t_tau - time);

        History history(rw->get_code(), rw->get_vertex(), time, rw->state_to_string(), "-", "Infected");
        rw->writeEvent(history);

        // increase number of contracted
        changeNumberContracted(num_Contracted + 1);
    } else
        v->sum_fail_encounters(infected);
}

void Simulator::process_walk(Event evt) {
    RandomWalk * rw = randomWalks.at(evt.randomwalk);

    // new time to walk
    double timeToWalk = rg.exponential(rw->get_lambda());
    rw->insertTimeWalking(timeToWalk);

    double t_walk = timeToWalk + time;

    // create new event and insert in events pq
    Event evt_new(t_walk, evt.randomwalk, EventType::Walk);
    events.push(evt_new);

    // get neighbor (vertex) to go
    int vIndex = rw->get_vertex();
    Vertex * srcVertex = vertices.at(vIndex);
    std::vector<unsigned int> neighbors = graph.get_vizinhos(vIndex + 1);
    int degree = graph.getGrau(vIndex + 1);
    int indexN = rg.uniform(0, degree); //uniform will not return b parameter
    int vertexNeighbor = neighbors.at(indexN);

    // delete rw in old vertex
    srcVertex->eraseRandomWalk(rw->get_rw_position());

    // vertex needs to be a pointer, if not gets a copy and don't change the real vertex in vector
    Vertex * dstVertex = vertices.at(vertexNeighbor - 1);
    randomWalks.at(evt.randomwalk)->set_rw_position(dstVertex->setRandomWalk(rw));
    randomWalks.at(evt.randomwalk)->set_vertex(vertexNeighbor - 1);

    //std::cout << "Vertex to go: " << vertexNeighbor -1 << std::endl;

    History history(rw->get_code(), rw->get_vertex(), time, rw->state_to_string(), eventToString(evt.type), "-");
    rw->writeEvent(history);


    if (rw->get_state() == State::Infected) {
        srcVertex->decreaseRwInfecteds(time);
        dstVertex->increaseRwInfecteds(time);
        infect(dstVertex, evt);
    } else if (rw->get_state() == State::Susceptible)
        beInfected(dstVertex, evt);
}

void Simulator::process_infect(Event evt) {
    // write snapshot of states
    writeNumberRwStatePerTime(evt.toString());

    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    rw->set_state(State::Infected);
    rw->setTimeStateChange(time, State::Infected);

    // generate time to recover
    double t_recover = time + rg.exponential(rw->get_gama());
    Event evt_new(t_recover, evt.randomwalk, EventType::Recover);
    events.push(evt_new);

    History history(rw->get_code(), rw->get_vertex(), time, rw->state_to_string(), eventToString(evt.type), "-");
    rw->writeEvent(history);
    rw->insertTimeInfected(t_recover - time);

    // increase number of vertices infected
    Vertex * vertex = vertices.at(rw->get_vertex());
    vertex->increaseRwInfecteds(time);

    change_infected_number(num_Infected + 1);
    num_Inf_Events--;

    // increase number of contracted
    changeNumberContracted(num_Contracted - 1);

    infect(vertex, evt);
}

void Simulator::process_heal(Event evt) {
    // write snapshot of states
    writeNumberRwStatePerTime(evt.toString());

    RandomWalk * rw = randomWalks.at(evt.randomwalk);
    rw->set_state(State::Susceptible);
    rw->setTimeStateChange(time, State::Susceptible);

    History history(rw->get_code(), rw->get_vertex(), time, rw->state_to_string(), eventToString(evt.type), "-");
    rw->writeEvent(history);

    change_infected_number(num_Infected - 1);

    Vertex * vertex = vertices.at(rw->get_vertex());
    vertex->decreaseRwInfecteds(time);

    beInfected(vertex, evt);
}

void Simulator::writeNumberRwStatePerTime(std::string evt) {
    // disable all file writing to be faster

    if (do_analysis) {
        std::ofstream arq;
        arq.open(fileNameNumberRandomWalkStates, std::ofstream::out | std::ofstream::app);
        arq << std::fixed << time << "," << num_Infected << "," << num_Inf_Events << "," << k - (num_Infected + num_Inf_Events) << "," << evt << std::endl;
        arq.close();
        
        infected_time.push_back(std::make_pair(time, num_Infected));
    }
}

void Simulator::write_infected_intervals() {
    if (do_analysis) {
        std::ofstream arq;
        arq.open(fileNameInfectInterval, std::ofstream::out | std::ofstream::app);

        for (int i = 0; i < k; i++)
            arq << i << "," << infected_intervals[i] / time << std::endl;

        arq.close();
    }
}

void Simulator::write_contracted_intervals() {
    if (do_analysis) {
        std::ofstream arq;
        arq.open(fileNameContractedInterval, std::ofstream::out | std::ofstream::app);

        for (int i = 0; i < k; i++)
            arq << i << "," << contracted_intervals[i] / time << std::endl;

        arq.close();
    }
}

void Simulator::write_susceptible_intervals() {

    if (do_analysis) {
        std::ofstream arq;
        arq.open(fileNameSusceptibleInterval, std::ofstream::out | std::ofstream::app);

        for (int i = 0; i < k; i++)
            arq << i << "," << susceptible_intervals[i] / time << std::endl;

        arq.close();
    }
}

void Simulator::write_infected_density() {
    Logger::Trace("Writing infected density over time");
    
    std::ofstream arq;
    arq.open(file_infected_density, std::ofstream::out | std::ofstream::app);

    for(auto it = infected_density_over_time.begin(); it != infected_density_over_time.end(); ++it)
        arq << (*it).first << ",             " << (*it).second << std::endl;

    arq.close();    
    
    Logger::Trace("End");
}

void Simulator::process() {

    //std::cout << "Begin epidemic" << std::endl;

    double old_percentual = 0.0;
    while (!events.empty()) {

        if (!continue_simulation_after_epidemic && (num_Inf_Events == 0 && num_Infected == 0))
            break;

        if (time > limit_time_epidemic)
            break;

        std::string evt_size = std::to_string(events.size());

        Event evt = get_top_event();
        
        // register fraction of infected
        if(evt.type == EventType::Infect)
            infected_density_over_time.push_back(std::make_pair(time, num_Infected/(double)k));
        
        time = evt.time;
        
        //Logger::Trace("Evento: " + evt.toString());
        //Logger::Trace("Events_0: " + evt_size);
        //Logger::Trace("Events_1: " + std::to_string(events.size()));  
        switch (evt.type) {
            case Walk:
                process_walk(evt);
                break;
            case Infect:
                process_infect(evt);
                break;
            case Recover:
                process_heal(evt);
                break;

            default:
                break;
        }
        //Logger::Trace("Events_2: " + std::to_string(events.size()));    
        if (do_analysis && evt.type != EventType::Walk)
            writeNumberRwStatePerTime(evt.toString());

        double percentual = time / limit_time_epidemic * 100;

        if ((int) percentual != (int) old_percentual) {
            Logger::Trace(std::to_string(percentual) + "%");
            //Logger::Trace(std::to_string(num_Infected) + ";" + std::to_string(num_Contracted) + ";" +  std::to_string(k - num_Contracted - num_Infected));
            //Logger::Trace("Events: " + std::to_string(events.size()));    
        }
        old_percentual = percentual;
    }

    //write_infected_density();

    //std::cout << "End epidemic" << std::endl;

    //std::cout << "#Tempo de simulação da epidemia: " << time << std::endl;

    // don't write to be fast, checking threshold
    if (do_analysis) {
//        for (int i = 0; i < randomWalks.size(); i++) {
//            randomWalks.at(i)->WriteResults();
//        }
//        
//        for (int i = 0; i < vertices.size(); i++) {
//            vertices.at(i)->writeTimeInfected(0.0);
//        }
//
//        write_infected_intervals();
//        write_contracted_intervals();
//        write_susceptible_intervals();
    }

}