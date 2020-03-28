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
	Fun��o do tipo void - Procedimento que faz o escalonamento da matriz.
	� executado da seguinte forma: 
	- recebe como par�metros uma matriz A[][], i e p: A � a matriz que ser� escalonada,
	  i � a linha a partir da qual ser�o feitas as transforma��es elementares, enquanto que
	  p � o �ndice para a coordenada (i, j) = (p, p) do piv�;
	- Chamadas - 1�: i=1 e j=0; 2�: i=2 e j=1...
*/
void triang1(double A[][n], int i, int p) //i: linhas que ser�o transformadas e p: piv� (A[p][p]
{
    double mi; //multiplicador
    int j = p;
    while (i < m){//transforma��es elementares de cada linha
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
    while (i < pi){ //faz as transforma��es elementares acima do piv�
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
    while (i < m){ //faz as transforma��es elementares abaixo do piv�
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

/*Fun��o do tipo int: retorna o �ndice i do maior elemento da linha, ou seja, do piv�.
Usada para aplicar o m�todo de Gauss sem pivota��o e com pivota��o parcial.*/
int busca_pivo(double A[][n], int p) //De in�cio, p � o �ndice para a linha e coluna do piv�. Nas itera��es, eventualmente p se tornar� o �ndice de alguma linha de baixo.
{
    int i, j = p; //j � o �ndice da coluna, que se mant�m fixo.
    double pivo = A[p][p]; //piv� inicial.

    for (i = p+1; i < m; i++){ //percurso linha a linha, come�ando pela imediatamente posterior � do piv�, em busca de um piv� maior
        if ( abs_mod(pivo) < abs_mod(A[i][j]) ){
            pivo = A[i][j]; // atualiza��o do valor do piv�, ou seja, encontrou um elemento maior que o piv� atual
            p = i; // atualiza��o do �ndice da linha do piv�
        }
    }
    return p;
}

/*Fun��o do tipo coord: retorna a coordenada (i, j) do maior elemento do sistema, isto �, do piv�.
Usada para aplicar o m�todo de Gauss com pivota��o total.*/
coord busca_Pivo(double A[][n], int pij)
{
    int i, j;
    coord p_ij; //para poder retornar um par ordenado (i, j) ao mesmo tempo
    p_ij.i = p_ij.j = pij; //coordenadas iniciais do piv�
    double pivo = A[pij][pij]; //piv� inicial
    for (i = pij; i < m; i++){ //percurso linha a linha
        for (j = 0; j < n-1; j++){ //percurso coluna a coluna. n-1, pois o piv� n�o pode ser um elemento da coluna de resposta
            if (abs_mod(pivo) < abs_mod(A[i][j])){
                pivo = A[i][j]; //atualiza��o do valor do piv�
                p_ij.i = i; //atualiza��o do �ndice i
                p_ij.j = j; //atualiza��o do �ndice j
            }
        }
    }
    return p_ij;
}

/*Fun��o do tipo double: retorna o valor absoluto (m�dulo) de um valor do tipo double.*/
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
    double x4, x3, x2, x1; // auxiliares para facilitar a manipula��o

    x4 = fixdig(A[3][4] / A[3][3], 4);
    X[3] = x4; //guardando

    x3 = fixdig(fixdig(A[2][4] - fixdig(x4*A[2][3], 4), 4) / A[2][2], 4);
    X[2] = x3;

    x2 = fixdig(fixdig(fixdig(A[1][4] - fixdig(x4*A[1][3], 4), 4) - fixdig(x3*A[1][2], 4), 4) / A[1][1], 4);
    X[1] = x2;

    x1 = fixdig(fixdig(fixdig(fixdig(A[0][4] - fixdig(x4*A[0][3], 4), 4) - fixdig(x3*A[0][2], 4), 4) - fixdig(x2*A[0][1], 4), 4) / A[0][0], 4);
    X[0] = x1;
}
