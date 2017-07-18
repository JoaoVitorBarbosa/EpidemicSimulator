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

    /// Generates a ring graph
    /// \param n
    /// \return 
    static ManipulaGrafoV Ring(int n) {
        ManipulaGrafoV graph;
        graph.inicializaGrafo(n);

        graph.incluirAresta(1, n, 1.0); // edge between zero and n.
        for (int i = 1; i < n; i++) {
            graph.incluirAresta(i, i + 1, 1.0); // always need to sum 1 (ManipulaGraph starts from 1)
        }

        return graph;
    }

    /// Generates a clique graph
    /// \param n vertices
    /// \return 
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

    /// Generates a bipartite graph
    /// \param n vertices of type 1
    /// \param n2 vertices of type 2
    /// \return 
    static ManipulaGrafoV Bipartite(int n, int n2) {
        ManipulaGrafoV graph;
        graph.inicializaGrafo(n + n2);
        for(int i=0; i < n; i++)
            for(int j= n; j < n+n2; j++)
                graph.incluirAresta(i+1, j+1, 1.0);
        
        return graph;
    }
};

#endif /* GRAPHGENERATOR_H */

