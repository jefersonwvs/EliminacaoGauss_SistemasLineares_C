#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "biblioteca_SistemasLineares.h"

double fixdig(double x, int dig)//arredonda
{
    int digInt, digFrac;
    if (x == 0.0)
        return 0.0;
    if (x < 0)
        digInt = floor((log10(-x))) + 1; // floor: piso
    else
        digInt = floor((log10(x))) + 1;
    digFrac = dig - digInt;
    return round(pow(10, digFrac) * x) / pow(10, digFrac); // round: arredondamento
}

/*	
	Função do tipo void - Procedimento que faz o escalonamento da matriz.
	É executado da seguinte forma: 
	- recebe como parâmetros uma matriz A[][], i e p: A é a matriz que será escalonada,
	  i é a linha a partir da qual serão feitas as transformações elementares, enquanto que
	  p é o índice para a coordenada (i, j) = (p, p) do pivô;
	- Chamadas - 1ª: i=1 e j=0; 2ª: i=2 e j=1...
*/
void triang1(double A[][n], int i, int p) //i: linhas que serão transformadas e p: pivô (A[p][p]
{
    double mi; //multiplicador
    int j = p;
    while (i < m){//transformações elementares de cada linha
        mi = fixdig((A[i][p]/A[p][p]), 4); //multiplicadores: A[1][0]/A[0][0]; A[2][1]/A[1][1]; ...
        printf("\n\tm%d = %lf/%lf ==> m%d = %lf", i+1, A[i][p], A[p][p], i+1, mi);
        A[i][p] = 0.0;
        for (j = p+1; j < n; j++){
            A[i][j] = fixdig(fixdig((-mi)*A[p][j], 4) + A[i][j], 4); //-multiplicador * linha pivotal + linha i
        }i++;
    }printf("\n");
}

void triang2(double A[][n], int a, coord pij)
{
    double mi;
    int i = a, j, pi = pij.i, pj = pij.j;
    while (i < pi){ //faz as transformações elementares acima do pivô
        mi = fixdig((A[i][pj]/A[pi][pj]), 4);
        printf("\n\tm%d = %lf/%lf ==> m%d = %lf", i+1, A[i][pj], A[pi][pj], i+1, mi);
        for (j = 0; j < n; j++){
            if (j == pj)
				A[i][j] = 0.0;
			else
				A[i][j] = fixdig(fixdig((-mi)*A[pi][j], 4) + A[i][j], 4);
        }
        i++;
    }
	
    i = pi+1;
    while (i < m){ //faz as transformações elementares abaixo do pivô
        mi = fixdig((A[i][pj]/A[pi][pj]), 4);
        printf("\n\tm%d = %lf/%lf ==> m%d = %lf", i+1, A[i][pj], A[pi][pj], i+1, mi);
        for (j = 0; j < n; j++){
            if (j == pj)
				A[i][j] = 0.0;
			else
				A[i][j] = fixdig(fixdig((-mi)*A[pi][j], 4) + A[i][j], 4);
        }
        i++;
    }
}

/*Função do tipo int: retorna o índice i do maior elemento da linha, ou seja, do pivô.
Usada para aplicar o método de Gauss sem pivotação e com pivotação parcial.*/
int busca_pivo(double A[][n], int p) //De início, p é o índice para a linha e coluna do pivô. Nas iterações, eventualmente p se tornará o índice de alguma linha de baixo.
{
    int i, j = p; //j é o índice da coluna, que se mantém fixo.
    double pivo = A[p][p]; //pivô inicial.

    for (i = p+1; i < m; i++){ //percurso linha a linha, começando pela imediatamente posterior à do pivô, em busca de um pivô maior
        if ( abs_mod(pivo) < abs_mod(A[i][j]) ){
            pivo = A[i][j]; // atualização do valor do pivô, ou seja, encontrou um elemento maior que o pivô atual
            p = i; // atualização do índice da linha do pivô
        }
    }
    return p;
}

/*Função do tipo coord: retorna a coordenada (i, j) do maior elemento do sistema, isto é, do pivô.
Usada para aplicar o método de Gauss com pivotação total.*/
coord busca_Pivo(double A[][n], int pij)
{
    int i, j;
    coord p_ij; //para poder retornar um par ordenado (i, j) ao mesmo tempo
    p_ij.i = p_ij.j = pij; //coordenadas iniciais do pivô
    double pivo = A[pij][pij]; //pivô inicial
    for (i = pij; i < m; i++){ //percurso linha a linha
        for (j = 0; j < n-1; j++){ //percurso coluna a coluna. n-1, pois o pivô não pode ser um elemento da coluna de resposta
            if (abs_mod(pivo) < abs_mod(A[i][j])){
                pivo = A[i][j]; //atualização do valor do pivô
                p_ij.i = i; //atualização do índice i
                p_ij.j = j; //atualização do índice j
            }
        }
    }
    return p_ij;
}

/*Função do tipo double: retorna o valor absoluto (módulo) de um valor do tipo double.*/
double abs_mod(double x)
{
    if (x == 0)
        return 0.0;
    else{
        if (x < 0)
            return -x;
        else
            return x;
    }
}

void troca_linha(double A[][n], int i, int j, int pi)
{
    double aux;
    while (j < n){
        aux = A[i][j];
        A[i][j] = A[pi][j];
        A[pi][j] = aux;
        j++;
    }
}

void soluc(double A[][n], double X[])
{
    double x4, x3, x2, x1; // auxiliares para facilitar a manipulação

    x4 = fixdig(A[3][4] / A[3][3], 4);
    X[3] = x4; //guardando

    x3 = fixdig(fixdig(A[2][4] - fixdig(x4*A[2][3], 4), 4) / A[2][2], 4);
    X[2] = x3;

    x2 = fixdig(fixdig(fixdig(A[1][4] - fixdig(x4*A[1][3], 4), 4) - fixdig(x3*A[1][2], 4), 4) / A[1][1], 4);
    X[1] = x2;

    x1 = fixdig(fixdig(fixdig(fixdig(A[0][4] - fixdig(x4*A[0][3], 4), 4) - fixdig(x3*A[0][2], 4), 4) - fixdig(x2*A[0][1], 4), 4) / A[0][0], 4);
    X[0] = x1;
}
