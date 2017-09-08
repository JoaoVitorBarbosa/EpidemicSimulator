/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Utils.cpp
 * Author: joao
 * 
 * Created on 8 de Setembro de 2017, 15:10
 */

#include "Utils.h"

Utils::Utils() {
}

Utils::Utils(const Utils& orig) {
}

Utils::~Utils() {
}

void Utils::createDirectory(std::string directory) {
    boost::filesystem::path dir(directory.c_str());
    if (boost::filesystem::create_directory(dir)) {
        std::cout << "Directory " << directory << " created successfully" << "\n";
    }
}


