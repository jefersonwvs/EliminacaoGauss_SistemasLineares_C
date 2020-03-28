#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <math.h>
#include "biblioteca_SistemasLineares.h"

int menu_principal();
int menu_metodos();
void inicializa_sistema(double[][n]);
void imprimir(double[][n], int, int, int*);
void inserir_elem(double[][n]);

int main()
{
    setlocale(LC_ALL, "Portuguese");

    double S[m][n], X[m]; //para um sistema qualquer; x1, x2, x3, x4
    int a, menu, p, linha;
    coord pij; //pij: coordenadas do pivô, para ser usado em pivotação total

    while (1){
        linha = 0;
        menu = menu_principal();
        if (menu == 1 || menu == 2){
            if (menu == 1){
                system("cls");
                printf("[1] Sistema nº 4:\n\n");
                inicializa_sistema(S); //copia S4 (sistema nº 4) para S
            }else{
                system("cls");
                printf("[2] Inserindo sistema:\n\n");
                inserir_elem(S);
            }
            imprimir(S, 0, 0, &linha);
            linha = 0;
            menu = menu_metodos();
            switch(menu){
                case 1:{ //sem pivotação
                    system("cls");
                    printf("Sistema:\n\n");
                    imprimir(S, 0, 0, &linha);
                    printf("[1] Aplicação do Método de Gauss sem pivotação\n\nEscalonando sistema...\n");
                    a = 0;
                    while (a < m-1){ //garante que zera somente abaixo das três primeiras colunas
                        triang1(S, a+1, a);
                        printf("\n");
                        imprimir(S, a+1, 0, &linha);
                        a++;
                    }
                    printf("...sistema escalonado:\n\n");
                    linha = 0;
                    imprimir(S, 0, 0, &linha);
                    soluc(S, X);
                    printf("Soluções:\t[ x1 = %lf | x2 = %lf | x3 = %lf | x4 = %lf ]\n\n", X[0], X[1], X[2], X[3]);
                }break;

                case 2:{ //com pivotação parcial;
                    system("cls");
                    printf("Sistema:\n\n");
                    imprimir(S, 0, 0, &linha);
                    printf("[2] Aplicação do Método de Gauss com pivotação parcial\n\nEscalonando sistema...\n");
                    a = 0;
                    while (a < m-1){
                        pij.i = busca_pivo(S, a);
                        pij.j = a;
                        printf("\n\tPivô: S[%d, %d] = %lf", pij.i, a, S[pij.i][a]);
                        triang2(S, a, pij);
                        //troca_linha(S, a, a, p);
                        //triang1(S, a+1, a);
                        troca_linha(S, a, 0, pij.i);
                        printf("\n");
                        imprimir(S, a+1, 0, &linha);
                        a++;
                    }
                    printf("...sistema escalonado:\n\n");
                    linha = 0;
                    imprimir(S, 0, 0, &linha);
                    soluc(S, X);
                    printf("Soluções:\t[ x1 = %lf | x2 = %lf | x3 = %lf | x4 = %lf ]\n\n", X[0], X[1], X[2], X[3]);
                }break;

                case 3:{ //com pivotação total;
                    system("cls");
                    printf("Sistema:\n\n");
                    imprimir(S, 0, 0, &linha);
                    printf("[3] Aplicação do Método de Gauss com pivotação total\n\nEscalonando sistema...\n");
                    a = 0;
                    while (a < m-1){
                        pij = busca_Pivo(S, a);
                        printf("\n\tPivô: S[%d, %d] = %lf", pij.i, pij.j, S[pij.i][pij.j]);
                        triang2(S, a, pij); printf("\n\n");
                        troca_linha(S, a, 0, pij.i);
                        imprimir(S, a+1, 0, &linha);
                        a++;
                    }
                    linha = 0;
                    printf("...sistema escalonado:\n\n");
                    linha = 0;
                    imprimir(S, 0, 0, &linha);
                }break;

                default: //retorna ao menu principal
                    system("cls");
                    goto atalho;
            }
            system("pause"); system("cls"); atalho:;
        }
        else{
            system("cls");
            printf("[3] ...encerrando programa.\n\n");
            break;
        }
    }
    return 3;
}

int menu_principal()
{
    int menu;
    printf("  ======Resolução de Sistemas Lineares====== \n"
           " | [1] Sistema nº 4;                        |\n"
           " | [2] Inserir sistema;                     |\n"
           " | [3] Encerrar programa.                   |\n"
           "  ========================================== \n >>> ");
    scanf("%d", &menu);
    while (menu < 1 || menu > 3){
        printf(" | [!] Opção inválida!                      |\n >>> "); scanf("%d", &menu);
    }
    printf("\n");
    return menu;
}

int menu_metodos()
{
    int menu;
    printf("  =============Método de Gauss============== \n"
           " | [1] Sem pivotação;                       |\n"
           " | [2] Com pivotação parcial;               |\n"
           " | [3] Com pivotação total;                 |\n"
           " | [4] Retornar ao menu principal.          |\n"
           "  ========================================== \n >>> ");
    scanf("%d", &menu);
    while (menu < 1 || menu > 4){
        printf(" | [!] Opção inválida!                      |\n >>> "); scanf("%d", &menu);
    }
    printf("\n");
    return menu;
}

void inicializa_sistema(double S[][n])
{
    /*double SQ[m][n] = {{0.065, -0.016, 0.013, 1.9, 9.575},          // Sistema 1 (Gelson)
                       {1.3, -0.31, 0.02, -1.9, -7.27},
                       {2.1, -0.26, -3.2, -0.49, 11.09},
                       {2.5, 0.011, 0.015, -1.2, -1.034}
                      };*/
    double SQ[m][n] = { {-0.082, -0.015, 0.0018, 1.4, -4.407},    // Sistema 4 (Jeferson Willian V. Silva)
                        {2.4, 0.42, 0.018, 0.94, 3.63},
                        {-1.4, -0.85, 2.5, 0.44, 8.68},
                        {-2.1, 0.0052, -0.0057, -0.92, -3.579} };
    /*double SQ[m][n] = {{-0.09, -1.47, 13.6, 18.6, -8.57},         // Sistema exemplo
                       {11.3, 2.56, -4.08, -5.67, -9.11},
                       {14.4, 13.1, 3.66, 15.9, 3.64},
                       {4.03, -3.89, 8.13, 0.850, 2.46}
                      };*/
    /*double SQ[m][n] = {{-0.093, 0.019, -0.016, -0.95, -1.18},     // Sistema 5 (Montese)
                       {-2.9, 0.59, -0.014, -1.6, -9.732},
                       {3.4, -1.1, 20, 1.5, -27.3},
                       {2.1, 0.0023, 0.013, 0.047, 4.212}
                      };*/
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            S[i][j] = SQ[i][j];
}

void inserir_elem(double A[][n])
{
    int i, j;
    for (i = 0; i < m; i++){
        for (j = 0; j < n; j++){
            printf("S[%d, %d] = ", i+1, j+1);
            scanf("%lf", &(A[i][j]));
        }
    }
}

void imprimir (double A[][n], int a, int b, int* l)
{
    int i, j;
    printf("\t\t (x1)\t\t (x2)\t\t (x3)\t\t (x4)\n");
    for (i = a; i < m; i++){
        printf("\t\t");
        for(j = b; j < n; j++){
            printf("[%lf]\t", A[i][j]);
        }
        (*l)++;
        printf("\t(L%d)\n", *l);
    }printf("\n");
}
