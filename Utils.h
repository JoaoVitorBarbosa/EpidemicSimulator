/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   Utils.h
 * Author: joao
 *
 * Created on 8 de Setembro de 2017, 15:10
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <stdio.h>
#include <boost/filesystem.hpp>
#include <string>

// create debug level
// create trace

class Utils {
private:
    static std::string pad(int num);
    
public:
    Utils();
    Utils(const Utils& orig);
    virtual ~Utils();
    
    /// Creates Directory to save files
    /// \param directory
    static void createDirectory(std::string directory);
    
    static std::string datetime_now_to_string();
    
    static std::string replace_dot_to_comma(double d);
    
private:

};

#endif /* UTILS_H */

