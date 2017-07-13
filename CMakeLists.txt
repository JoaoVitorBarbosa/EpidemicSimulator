cmake_minimum_required(VERSION 3.7)
project(Epidemic_Simulator)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        nbproject/private/c_standard_headers_indexer.c
        nbproject/private/cpp_standard_headers_indexer.cpp
        EventProcessor.h
        main.cpp
        ManipulaGrafo.cpp
        ManipulaGrafo.h
        RandomGenerator.cpp
        RandomGenerator.h
        RandomWalk.cpp
        RandomWalk.h
        Simulator.cpp
        Simulator.h
        Vertex.cpp
        Vertex.h)

add_executable(Epidemic_Simulator ${SOURCE_FILES})