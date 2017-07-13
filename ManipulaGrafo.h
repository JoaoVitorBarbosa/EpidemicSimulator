#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <stdlib.h>


#ifndef MANIPULAGRAFO_H
#define MANIPULAGRAFO_H

// Classe Base

class ManipulaGrafo{

public:

	//***********************************************************
	//***************** QUESTAO 1 - ENTRADA *********************
	//***********************************************************

	ManipulaGrafo();									//done
	~ManipulaGrafo();									//done

	//***********************************************************
	//***************** QUESTAO 2 - SAIDA ***********************
	//***********************************************************

	void gerarTxt();		// Saida em *.txt					//done

	//***********************************************************
	//***************** QUESTAO 3 - REPRESENTACAO ***************
	//***********************************************************

	//vector<int> getVetor(unsigned int); // retorna o vetor de adjacencias de vertice desejado

	//***********************************************************
	//***************** QUESTAO 4 - ALGORITMOS DE BUSCA *********
	//***********************************************************

	std::vector<std::vector<unsigned int>> BFS(unsigned int);	//done
	std::vector<std::vector<unsigned int>> DFS(unsigned int);	//done
    void ArvoreBFS(unsigned int);								//done
	void ArvoreDFS(unsigned int);								//done

	//***********************************************************
	//***************** QUESTAO 5 - COMPONENTES CONEXOS *********
	//***********************************************************

	void getCompConexo();										//done

	//***********************************************************
	//***************** OUTRAS FUNCIONALIDADES ******************
	//***********************************************************

	void get_diametro();
	unsigned int get_num_vertices(); // retorna o numero de vertices no grafo
	unsigned int get_num_arestas(); // retorna o numero de arestas no grafo
	std::string get_docTexto();	

        int getGrau(int index);
	virtual std::vector<unsigned int> get_vizinhos(unsigned int) = 0; //função que retorna um vetor de vizinhos de determinado vértice

        unsigned int num_vertices; //quantidade de vertices
	unsigned long int num_arestas; // quantidade de arestas
	std::vector<unsigned int> graus; // guarda os graus de todos os vertices
        void lerArquivo(std::string, bool); // le o arquivo de texto e constroi um grafo a partir dele //done
        void inicializaGrafo(int n);
        
        unsigned int grau(unsigned int); //retorna o grau de determinado vértice
	virtual void preparaEstrutura() = 0; // prepara a estrutura de dados (vetor ou matriz) para receber valores

	virtual void incluirAresta(unsigned int, unsigned int, double = 1) = 0; // inclui uma aresta no grafo

	std::string docTexto; // documento de texto usado para construir o grafo
};

#endif


//Classe herdeira para matrizes

#ifndef MANIPULAGRAFOM_H
#define MANIPULAGRAFOM_H

class ManipulaGrafoM : public ManipulaGrafo{
public:
	ManipulaGrafoM(std::string arquivo, bool compeso) : ManipulaGrafo(){
		docTexto = arquivo;
		lerArquivo(arquivo, compeso);
	};

	~ManipulaGrafoM(){};

	std::vector<unsigned int> get_vizinhos(unsigned int); //função que retorna um vetor de vizinhos de determinado vértice

	void preparaEstrutura();

	void incluirAresta(unsigned int, unsigned int, double = 1); // inclui uma aresta no grafo

	std::vector<std::vector<bool>> matriz;
};

#endif;


//Classe herdeira para vetores de adjacências

#ifndef MANIPULAGRAFOV_H
#define MANIPULAGRAFOV_H

class ManipulaGrafoV : public ManipulaGrafo{
public:
    	ManipulaGrafoV(){ 
	};
    
	ManipulaGrafoV(std::string arquivo, bool compeso) : ManipulaGrafo(){ 
		docTexto = arquivo;
		lerArquivo(arquivo, compeso);
	};

	~ManipulaGrafoV(){};

	std::vector<unsigned int> get_vizinhos(unsigned int); //função que retorna um vetor de vizinhos de determinado vértice

	void preparaEstrutura();

	void incluirAresta(unsigned int, unsigned int, double = 1); // inclui uma aresta no grafo

	std::vector<std::vector<unsigned int>> vetorAdj;
};

#endif;