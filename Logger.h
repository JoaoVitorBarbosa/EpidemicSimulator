/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Debug.h
 * Author: joao
 *
 * Created on 2 de Janeiro de 2018, 13:46
 */

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <string>

enum LogLevel {
    Error = 0,
    Warn,
    Info,
    Debug,
    Trace
};

class Logger {
private:

    static void Write(LogLevel level, std::string str) {
        std::cout << level_to_string(level) << ": " << str << std::endl;
    }

    static std::string level_to_string(LogLevel level) {
        switch (level) {
            case LogLevel::Debug:
                return "DEBUG";
            case LogLevel::Error:
                return "ERROR";
            case LogLevel::Info:
                return "INFO";
            case LogLevel::Trace:
                return "TRACE";
            case LogLevel::Warn:
                return "WARN";
        }
    }

public:
    static LogLevel level;

    Logger() {
        level = LogLevel::Error;
    }

    Logger(LogLevel _level) {
        level = _level;
    }

    static void Error(std::string str) {
        std::cout << level_to_string(LogLevel::Error) << ": " << str << std::endl;
    }

    static void Warn(std::string str) {
        if (level != LogLevel::Error)
            Write(LogLevel::Warn, str);
    }

    static void Info(std::string str) {
        if (level != LogLevel::Warn && level != LogLevel::Error)
            Write(LogLevel::Info, str);
    }

    static void Debug(std::string str) {
        if (level == LogLevel::Debug)
            Write(LogLevel::Debug, str);
    }

    static void Trace(std::string str) {
        if (level == LogLevel::Trace)
            Write(LogLevel::Trace, str);
    }
};

#endif /* DEBUG_H */

