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

std::string Utils::datetime_now_to_string() {
    time_t now = std::time(0);
    tm *ltm = localtime(&now);
    std::string time_str = Utils::pad(1900 + ltm->tm_year) + Utils::pad(1 + ltm->tm_mon) + Utils::pad(ltm->tm_mday) + Utils::pad(ltm->tm_hour) + Utils::pad(ltm->tm_min) + Utils::pad(1 + ltm->tm_sec);
    return time_str;
}

void Utils::createDirectory(std::string directory) {
    boost::filesystem::path dir(directory.c_str());
    if (boost::filesystem::create_directory(dir)) {
        //std::cout << "Directory " << directory << " created successfully" << "\n";
    }
}

std::string Utils::pad(int num)
{
    std::string aux = num < 10 ? "0" + std::to_string(num) : std::to_string(num);
    return aux;
}