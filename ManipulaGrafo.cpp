#include "ManipulaGrafo.h"

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ************************************************************
 ************************************************************
 ********************* ManipulaGrafo ************************
 ************************************************************
 ************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */



ManipulaGrafo::ManipulaGrafo() {
    num_vertices = 0;
    num_arestas = 0;
}

ManipulaGrafo::~ManipulaGrafo() {
}


//***********************************************************
//***************** QUESTAO 2 - SAIDA ***********************
//***********************************************************

void ManipulaGrafo::gerarTxt() {
    // Saida em *.txt 

    unsigned int n = get_num_vertices();
    unsigned int m = get_num_arestas();

    double d_medio = 2 * (double) m / (double) n;

    std::ofstream arq;
    const char* arquivo = "Saida_do_Programa.txt";
    arq.open(arquivo);
    arq << "# n = " << n
            << "\n# m = " << m
            << "\n# d_medio = " << d_medio << "\n";


    unsigned int g_max = 0; // grau maximo do vertice

    // descobre o grau maximo do vertice
    for (unsigned int i = 1; i <= n; i++) {

        if (grau(i) > g_max)
            g_max = grau(i);
    }

    std::vector<unsigned int> distribuicaoEmpirica(g_max + 1, 0); // guarda o numero de ocorrencias de cada grau de vertice do grafo
    // cada elemento eh o numero de ocorrencias de cada grau
    // cada grau eh o indice no vetor que o representa.
    //Ex.: grau 0 = distribuicaoEmpirica[0]


    // conta a quantidade de cada grau possivel aos vertices 

    for (unsigned int i = 1; i <= n; i++) {

        distribuicaoEmpirica.at(grau(i))++;

    }

    // escreve a distribuicao empirica no arquivo

    double distribuicao = 0.0; // variavel usada para calcular a distribuicao empirica

    for (unsigned int i = 0; i <= g_max; i++) {

        distribuicao = (double) distribuicaoEmpirica.at(i) / (double) m;

        arq << i << " " << distribuicao << "\n";
    }

    arq.close();
}

int ManipulaGrafo::getGrau(int index) {
    return graus.at(index - 1);
}

//***********************************************************
//***************** QUESTAO 3 - REPRESENTACAO ***************
//***********************************************************


//***********************************************************
//***************** QUESTAO 4 - ALGORITMOS DE BUSCA *********
//***********************************************************

std::vector<std::vector<unsigned int>> ManipulaGrafo::BFS(unsigned int v) {

    std::vector<std::vector<unsigned int>> ag; //vetor que conterá a arvore geradora. Usado pela função ArvoreBFS(unsigned int) 
    //para escrever o arquivo de texto
    std::vector<bool> vert(get_num_vertices(), false); // cria vetor contendo todos os vértices, desmarcados.
    // cria a fila
    std::queue<std::vector<unsigned int>> Q;

    std::vector<unsigned int> paiV; // vetor que mostra o vertice 'v', seu pai na arvore  e seu nível
    // vetores parecidos sao usados no loop for

    paiV.resize(3);
    paiV.at(0) = v; // vetor(0) - id do vertice
    paiV.at(1) = 0; // vetor(1) - pai do vertice
    paiV.at(2) = 0; // vetor(2) - nivel do vertice

    Q.push(paiV);

    vert.at(v - 1) = true;

    ag.push_back(paiV);

    while (!Q.empty()) {

        std::vector<unsigned int> u = Q.front(); //pegar primeiro elemento
        Q.pop();

        // checar todas as arestas dizendo que o vertice em questao eh pai de todos estes (loop for) e têm nivel = nivel + 1
        // checar se ja estao marcados antes de colocar na fila

        std::vector<unsigned int> vizinhos_u(grau(u[0]), 0); // define vetor que guardará todos os vizinhos do vertice u

        vizinhos_u = get_vizinhos(u[0]); // guarda todos os vizinhos do vertice u // <---------------------------------- o problema esta aqui!

        for (unsigned int k = 0; k < vizinhos_u.size(); k++) {
            if (vert[vizinhos_u[k] - 1] == false) {
                vert[vizinhos_u[k] - 1] = true; // marca o vertice

                std::vector<unsigned int> w(3, 0); //prepara o vertice para ser incluido na fila
                w.resize(3);
                w[0] = vizinhos_u[k]; //id
                w[1] = u[0]; //pai
                w[2] = u[2] + 1; //nivel

                Q.push(w); //inclui vizinho de u na fila

                ag.push_back(w); // inclui vizinho de u na arvore geradora

            }
        }
    }
    return ag;
}

std::vector<std::vector<unsigned int>> ManipulaGrafo::DFS(unsigned int v) {

    std::vector<std::vector<unsigned int>> ag; //vetor que conterá a arvore geradora. Usado pela função ArvoreBFS(unsigned int) 
    //para escrever o arquivo de texto

    std::vector<bool> vert(get_num_vertices(), false); // cria vetor contendo todos os vértices, desmarcados.
    // cria a fila
    std::stack<std::vector<unsigned int>> P;

    std::vector<unsigned int> paiV; // vetor que mostra o vertice 'v' e seu pai na arvore = 0
    // vetores parecidos sao usados no loop for

    paiV.resize(3);
    paiV.at(0) = v; // vetor(0) - id do vertice
    paiV.at(1) = 0; // vetor(1) - pai do vertice
    paiV.at(2) = 0; // vetor(2) - nivel do vertice

    P.push(paiV);

    while (!P.empty()) {
        std::vector<unsigned int> u = P.top(); //pegar primeiro elemento
        P.pop(); // remover o elemento da pilha

        if (vert[u[0] - 1] == false) { //se o vertice não estiver marcado
            /*Relembrando: A estrutura dos vetores eh
            u[0] - id do vertice
            u[1] - pai do vertice
            u[2] - nivel do vertice
             */

            vert.at(u.at(0) - 1) = true; // marca o vertice

            ag.push_back(u);

            std::vector<unsigned int> vizinhos_u(grau(u[0]), 0); // define vetor que guardará todos os vizinhos do vertice u

            vizinhos_u = get_vizinhos(u[0]); // guarda todos os vizinhos do vertice u

            for (unsigned int j = 0; j < vizinhos_u.size(); j++) {

                std::vector<unsigned int> w; // novo vertice a ser incluido na pilha
                w.resize(3);
                w.at(0) = vizinhos_u[j]; // id
                w.at(1) = u.at(0); // pai
                w.at(2) = u.at(2) + 1; // nivel

                P.push(w); // inclui adjacencia na pilha sem marca-la

                ag.push_back(w);
            }
        }
    }
    return ag;
}

void ManipulaGrafo::ArvoreBFS(unsigned int v) {

    std::ofstream ag("Arvore_Geradora_BFS.txt", std::ios::out); // arquivo onde será escrita a arvore geradora

    std::vector<std::vector<unsigned int>> arvore = BFS(v);

    // arvore[x][0] == v
    // arvore[x][1] == pai de v
    // arvore[x][2] == nivel de v

    ag << "\n\nArvore geradora da BFS: \n\n" << v << ":				nivel = " << arvore[0][2]
            << "	Raiz da Árvore \n";

    for (unsigned int i = 1; i < arvore.size(); i++) {

        ag << arvore[i][0] << ":				nivel = " << arvore[i][2]
                << "	pai: " << arvore[i][1] << " \n";

    }
    ag.close();
}

void ManipulaGrafo::ArvoreDFS(unsigned int v) {

    //unsigned int nivel = 0; // nivel do vertice na arvore
    std::ofstream ag("Arvore_Geradora_DFS.txt", std::ios::out); // arquivo onde a arvore sera descrita

    std::vector<std::vector<unsigned int>> arvore = DFS(v);

    ag << "\n\nArvore geradora da DFS: \n\n";

    ag << arvore[0][0] << ":				nivel = " << arvore[0][2]
            << "	pai: Não possui \n";

    for (unsigned int k = 1; k < arvore.size(); k++) {

        ag << arvore[k][0] << ":				nivel = " << arvore[k][2]
                << "	pai: " << arvore[k][1] << "\n";

    }
    ag.close();
}


/*// <função private>

unsigned int ManipulaGrafo::ModBFS(unsigned int i) // função private
{	
        unsigned int diametro = 0; // Diametro do grafo

        std::vector<unsigned int> vert(this->num_vertices, false);
        std::vector<unsigned int> nivel(this->num_vertices, 0);
        std::vector<unsigned int> pai(this->num_vertices, 0);

        std::queue<unsigned int> fila;

        vert.at(i - 1) = true;// vertice marcado
        nivel.at(i - 1) = 0;
        pai.at(i - 1) = 0;

        fila.push(i);

        while (!fila.empty())
        {
                int v = fila.front();
                fila.pop();

                for (unsigned int j = 0; j < arestas.at(v - 1).adjacencia.size(); j++)
                        {

                                if (!vert.at(arestas.at(v - 1).adjacencia.at(j) - 1))
                                {
                                        // arestas.at(v - 1).adjacencia.at(j) é o vertice que a BFS encontrou

                                        vert.at(arestas.at(v - 1).adjacencia.at(j) - 1) = true; //vertice visitado
                                        fila.push(arestas.at(v - 1).adjacencia.at(j)); // inclui vertice na fila
                                        pai.at(arestas.at(v - 1).adjacencia.at(j) - 1) = v; // diz quem é o pai do vertice
                                        nivel.at(arestas.at(v - 1).adjacencia.at(j) - 1) = nivel.at(v-1) + 1; // diz qual é o nivel do vertice

                                        unsigned int diamet = nivel.at(arestas.at(v - 1).adjacencia.at(j) - 1);

                                        if (diametro < diamet)
                                                diametro = diamet;
                                }
                        }
        }


        return diametro;

        }*/
// </função private>
//ModBFS()


//***********************************************************
//***************** QUESTAO 5 - COMPONENTES CONEXOS *********
//***********************************************************

void ManipulaGrafo::getCompConexo() {

    std::ofstream ag("Componentes_Conexos.txt", std::ios::out);
    unsigned int count = 1; // contador de componentes conexos;

    std::vector<bool> vert(get_num_vertices(), false); // cria vetor contendo todos os vértices, desmarcados.

    for (unsigned int i = 0; i < get_num_vertices(); i++) {
        if (vert[i] == false) {

            std::vector<std::vector<unsigned int>> arvore = BFS(i + 1);

            ag << "\nComponente Conexo " << count << ":\n\n";

            for (unsigned int j = 0; j < arvore.size(); j++) {
                vert[arvore[j][0] - 1] = true; //marca os vertices ja explorados
                ag << arvore[j][0] << ":				nivel = " << arvore[j][2]
                        << "	pai: " << arvore[j][1] << " \n";

            }
            count++;
        }
    }
    ag.close();
}


//***********************************************************
//***************** OUTRAS FUNCIONALIDADES ******************
//***********************************************************

void ManipulaGrafo::get_diametro() {

    std::vector<std::vector<unsigned int>> arvore; // arvore da BFS

    unsigned int diam = 0; // diametro de cada arvore
    unsigned int diametro = 0; // diametro do grafo

#pragma omp parallel for schedule(dynamic)    //Codigo simples para paralelizar, o schedule(dynamic) é para quando o tempo de resposta de cada tarefa pode variar.
    for (int i = 1; i <= get_num_vertices(); i++) {
        arvore = BFS(i);
        unsigned int size = arvore.size();
        diam = arvore[size - 1][2]; // diametro da BFS = nível do ultimo vértice na arvore
        //diam = arvore.back()[2];
#pragma omp critical(dataupdate)   //nesse bloco, somente um thread por vez pode acessar. Assim quando cada thread acaba seu pacote de BFS, faz a comparação. Cada thread terá seu diametro maximo, e tiro o diametro maximo dos threads
        {
            if (diametro < diam)
                diametro = diam;
        }
    }
    std::cout << "\nDiametro = " << diametro << "\n";
}

unsigned int ManipulaGrafo::get_num_vertices() {
    return this->num_vertices;
}

unsigned int ManipulaGrafo::get_num_arestas() {
    return this->num_arestas;
}

unsigned int ManipulaGrafo::grau(unsigned int v) {
    return graus.at(v - 1);
}

std::string ManipulaGrafo::get_docTexto() {
    return docTexto;
}

//private:

// deve ser usado quando não criar grafo a partir de arquivos

void ManipulaGrafo::inicializaGrafo(int n) {
    this->num_vertices = n;
    graus.assign(this->num_vertices, 0); // muda o tamanho do vetor graus para o numero de vertices e inicia todas as posições com 0
    preparaEstrutura(); // prepara a estrutura de dados para receber as arestas   
}

void ManipulaGrafo::lerArquivo(std::string arquivo, bool compeso) {

    //Esta função lê o arquivo linha a linha, 
    //selecionando os vértices e o peso da aresta corresponente.

    //Após isto, ela adiciona as arestas ao grafo utilizando
    //o método virtual incluirAresta(unsigned int, unsigned int, double)
    //que deve ser explicitado nas classes herdeiras.

    std::string linha; // linha que será lida do arquivo.
    std::string espaco = " "; // usado para contornar o problema de "char != const char"
    unsigned int count = 0; // contador de linhas

    std::ifstream arq;

    arq.open(arquivo/*, arq.in*/); // abre o arquivo


    while (!arq.eof()) {
        // enquanto o arquivo não acabar

        getline(arq, linha); // pegar uma linha
        count++;

        if (count == 1) {

            this->num_vertices = atof(linha.c_str()); // converte a string linha em inteiro
            graus.assign(this->num_vertices, 0); // muda o tamanho do vetor graus para o numero de vertices e inicia todas as posições com 0
            preparaEstrutura(); // prepara a estrutura de dados para receber as arestas
        }

        /*Agora, o programa lerá o arquivo para incluir as arestas*/

        // estes sao os "vertices" usados para incluir no grafo e na aresta. Devem ser convertidos para int
        std::string str_vertice1;
        unsigned int vertice1 = 0;
        std::string str_vertice2;
        unsigned int vertice2 = 0;
        std::string str_peso;
        double peso = 1.0; // Peso padrão é igual a 1. Isto faz com que a função funcione mesmo que o grafo não possua pesos em suas arestas


        if (count > 1) {

            // **** PESO ****
            if (compeso) {
                int firstSpace = linha.find_first_of(' ');
                if (firstSpace == -1)
                    break;

                int lastSpace = linha.find_last_of(' ');

                std::string vertex1 = linha.substr(0, firstSpace);
                std::string vertex2 = linha.substr(firstSpace++, lastSpace);
                str_peso = linha.substr(lastSpace++, linha.length());

                vertice1 = atof(vertex1.c_str());
                vertice2 = atof(vertex2.c_str());
                peso = atof(str_peso.c_str()); // converte a string str_peso em double

                // **** FIM PESO ****
            } else {

                int spaceIndex = linha.find(' ');
                if (spaceIndex == -1)
                    break;

                std::string vertex1 = linha.substr(0, spaceIndex);
                std::string vertex2 = linha.substr(spaceIndex++, linha.length());

                vertice1 = atof(vertex1.c_str());
                vertice2 = atof(vertex2.c_str());
            }

            incluirAresta(vertice1 + 1, vertice2 + 1, peso); // inclui a aresta	
        }
    }
    arq.close();
}

void ManipulaGrafo::incluirAresta(unsigned int v1, unsigned int v2, double peso) {
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ************************************************************
 ************************************************************
 ********************* ManipulaGrafoM ***********************
 ************************************************************
 ************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

//private:

void ManipulaGrafoM::preparaEstrutura() {
    std::vector<bool> vizinhos(get_num_vertices(), false);
    matriz.assign(get_num_vertices(), vizinhos);
}

std::vector<unsigned int> ManipulaGrafoM::get_vizinhos(unsigned int v) {

    std::vector<unsigned int> vizinhos;

    for (unsigned int i = 0; i < get_num_vertices(); i++) {
        if (matriz[v - 1][i] == true) {
            vizinhos.push_back(i + 1);
        }
    }
    return vizinhos;
}

void ManipulaGrafoM::incluirAresta(unsigned int v1, unsigned int v2, double peso) {

    if (!matriz[v1 - 1][v2 - 1] && !matriz[v2 - 1][v1 - 1]) {

        matriz[v1 - 1][v2 - 1] = peso; // torna v1 vizinho de v2
        matriz[v2 - 1][v1 - 1] = peso; // torna v2 vizinho de v1

        graus[v1 - 1]++; //incrementa o grau de v1
        graus[v2 - 1]++; //incrementa o grau de v2

        num_arestas++;
    }

}





/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ************************************************************
 ************************************************************
 ********************* ManipulaGrafoV ***********************
 ************************************************************
 ************************************************************
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

//private:

void ManipulaGrafoV::preparaEstrutura() {
    vetorAdj.resize(num_vertices);
}

std::vector<unsigned int> ManipulaGrafoV::get_vizinhos(unsigned int v) {
    return vetorAdj[v - 1];
}

void ManipulaGrafoV::incluirAresta(unsigned int v1, unsigned int v2, double peso) {

    bool permissao = peso; //permissao para adicionar a aresta

    //checa se v1 e v2 já são vizinhos
    for (unsigned int i = 0; i < vetorAdj[v1 - 1].size(); i++) {
        if (vetorAdj[v1 - 1][i] == v2) {
            permissao = false;
        }
    }

    //inclui a aresta
    if (permissao) {
        vetorAdj[v1 - 1].push_back(v2); // torna v1 vizinho de v2
        vetorAdj[v2 - 1].push_back(v1); // torna v2 vizinho de v1

        graus[v1 - 1]++; //incrementa o grau de v1
        graus[v2 - 1]++; //incrementa o grau de v2

        num_arestas++;
    }

}