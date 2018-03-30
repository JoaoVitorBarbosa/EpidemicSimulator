/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Event.h
 * Author: joao
 *
 * Created on 30 de Maio de 2017, 23:30
 */

#ifndef EVENT_H
#define EVENT_H

#include<iostream>
#include<string>

/// Represents type of events stored in Event list
enum EventType {
    Walk,
    Recover,
    Infect
};

class Event {
public:

    // time to event occur
    double time;
    // rw index in vector
    int randomwalk;
    EventType type;

    std::string toString() {
        switch (type) {
            case EventType::Infect:
                return "Infect";
            case EventType::Recover:
                return "Recover";
            case EventType::Walk:
                return "Walk";
        }

    }

    Event() {
    }

    Event(double _time, int _index, EventType _type) : time(_time), randomwalk(_index), type(_type) {
    }

    
    /// Overload operator to use Event in a min priority queue.
    friend bool operator<(const Event& e, const Event& r) {
        return e.time
                > r.time; // operator > inverted to change order
    }
};

#endif /* EVENT_H */


