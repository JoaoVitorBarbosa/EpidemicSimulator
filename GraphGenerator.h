/*
 * Project develop by João Vitor Barbosa Tavares
 *   * 
 */

/* 
 * File:   GraphGenerator.h
 * Author: joao
 *
 * Created on 11 de Julho de 2017, 14:35
 */

#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include "ManipulaGrafo.h"

class GraphGenerator {
public:
    
    static ManipulaGrafoV Ring(int n) {
        ManipulaGrafoV graph;
        graph.inicializaGrafo(n);
        
        graph.incluirAresta(1, n, 1.0); // edge between zero and n.
        for (int i = 1; i < n; i++)
        {
            graph.incluirAresta(i, i + 1, 1.0); // always need to sum 1 (ManipulaGraph starts from 1)
        }
        
        return graph;
    }

    static ManipulaGrafoV Clique(int n) {

        ManipulaGrafoV graph;
        graph.inicializaGrafo(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j)
                    graph.incluirAresta(i + 1, j + 1, 1.0);

        return graph;
    }

    static void Torus(int n) {

    }
};

#endif /* GRAPHGENERATOR_H */

