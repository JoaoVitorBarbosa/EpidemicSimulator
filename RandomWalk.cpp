/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

#include "RandomWalk.h"
#include "Simulator.h"

RandomWalk::RandomWalk() {
}

RandomWalk::RandomWalk(int _vertex, double _lambda, double _gama, double _tau, State _state, int _code, std::string _outputDir) {
    do_analysis = true;
    vertex = _vertex;
    lambda = _lambda;
    gama = _gama;
    tau = _tau;
    state = _state;
    code = _code;
    output_dir = _outputDir;

    std::string prefix = "RW" + std::to_string(code);
    file_name_history = output_dir + "/" + prefix + "_change  .txt";
    file_name_sample_results = output_dir + "/" + prefix + "_Results.txt";
    file_name_walking_times = output_dir + "/" + prefix + "_walkingTimes.txt";

    // commented
    if (do_analysis) {
        std::ofstream arq;

        arq.open(file_name_history);
        arq << "# Random Walk: " << code << std::endl;
        arq << "# Parameters: Lambda: " << lambda << ", Gama: " << gama << ", Tau: " << tau << '\n';
        arq << "# time,event,vertex,state,effect,infected_by" << '\n';

        arq.close();

        arq.open(file_name_sample_results);
        arq << "# Parameters: Lambda: " << lambda << ", Gama: " << gama << ", Tau: " << tau << '\n';
        arq << "# time,state_changed" << '\n';
        arq.close();

        arq.open(file_name_walking_times);
        arq << "# Walking Times" << std::endl;
        arq.close();
    }
}

void RandomWalk::set_vertex(int v) {
    vertex = v;
}

int RandomWalk::get_vertex() {
    return vertex;
}

void RandomWalk::set_code(int i) {
    code = i;
}

int RandomWalk::get_code() {
    return code;
}

void RandomWalk::set_state(State _state) {
    state = _state;
}

State RandomWalk::get_state() {
    return state;
}

void RandomWalk::set_lambda(double _lambda) {
    lambda = _lambda;
}

double RandomWalk::get_lambda() {
    return lambda;
}

void RandomWalk::set_gama(double _gama) {
    gama = _gama;
}

double RandomWalk::get_gama() {
    return gama;
}

void RandomWalk::set_tau(double _tau) {
    tau = _tau;
}

double RandomWalk::get_tau() {
    return tau;
}

void RandomWalk::set_rw_position(std::list<RandomWalk*>::iterator it) {
    rw_iterator = it;
}

std::list<RandomWalk*>::iterator RandomWalk::get_rw_position() {
    return rw_iterator;
}

bool RandomWalk::is_infected() {
    return state == Contracted;
}

std::string RandomWalk::state_to_string() {
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
    // commented
    if (do_analysis) {
        //std::cout << "RW: " << code << " ; Vertex: " << vertex << " ; Time: " << time << " ; Event: " << history.event << " ; State: " << state_to_string() << std::endl;
        writeFile(history);
    }
}

void RandomWalk::writeFile(History history) {
    // commented
    if (do_analysis) {
        std::ofstream arq;
        arq.open(file_name_history, std::ofstream::out | std::ofstream::app);
        arq << std::fixed << history.time << "," << history.event << "," << history.vertex << "," << history.state << "," << history.effect << "," << history.get_num_rw_infected() << '\n';
        arq.close();

        if (history.event == "Walk") {
            arq.open(file_name_walking_times, std::ofstream::out | std::ofstream::app);
            arq << std::fixed << history.time << std::endl;
            arq.close();
        }
    }
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
    // commented

    if (do_analysis) {
        std::ofstream arq;
        arq.open(file_name_sample_results, std::ofstream::out | std::ofstream::app);

        arq << std::fixed << "# Dados Infectado - Média: " << boost::accumulators::mean(intervalsInfected)
                << " Mediana: " << boost::accumulators::median(intervalsInfected)
                << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsInfected)) << '\n';

        arq << std::fixed << "# Dados Contraído - Média: " << boost::accumulators::mean(intervalsContracted)
                << " Mediana: " << boost::accumulators::median(intervalsContracted)
                << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsContracted)) << '\n';

        arq << std::fixed << "# Dados Caminhar - Média: " << boost::accumulators::mean(intervalsWalking)
                << " Mediana: " << boost::accumulators::median(intervalsWalking)
                << " Desvio Padrão: " << std::sqrt(boost::accumulators::variance(intervalsWalking)) << '\n';

        arq.close();
    }
}

void RandomWalk::setTimeStateChange(double _t, int _s) {
    // commented

    if (do_analysis) {
        state_timestamps.push_back(std::make_pair(_t, _s));

        std::ofstream arq;
        arq.open(file_name_sample_results, std::ofstream::out | std::ofstream::app);

        arq << std::fixed << _t << "," << _s << '\n';

        arq.close();
    }
}

void RandomWalk::drawStates() {
    Gnuplot gp;

    std::vector<std::vector<std::pair<double, int> > > suscetibles_segments;
    std::vector<std::vector<std::pair<double, int> > > contracteds_segments;
    std::vector<std::vector<std::pair<double, int> > > infecteds_segments;

    for (int i = 1; i < state_timestamps.size(); i++) {
        // time infected
        if (state_timestamps[i].second == State::Susceptible) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(state_timestamps[i - 1].first, 3));
            segment.push_back(std::make_pair(state_timestamps[i].first, 3));
            infecteds_segments.push_back(segment);
        }

        // time suscetible
        if (state_timestamps[i].second == State::Contracted) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(state_timestamps[i - 1].first, 1));
            segment.push_back(std::make_pair(state_timestamps[i].first, 1));
            suscetibles_segments.push_back(segment);
        }

        // time contracted
        if (state_timestamps[i].second == State::Infected) {
            std::vector<std::pair<double, int> > segment;
            segment.push_back(std::make_pair(state_timestamps[i - 1].first, 2));
            segment.push_back(std::make_pair(state_timestamps[i].first, 2));
            contracteds_segments.push_back(segment);
        }
    }

    int max = state_timestamps[state_timestamps.size() - 1].first;

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