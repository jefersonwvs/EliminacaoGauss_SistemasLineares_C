#ifndef BIBLIOTECA_SISTEMASLINEARES_H_INCLUDED
#define BIBLIOTECA_SISTEMASLINEARES_H_INCLUDED

#define m 4
#define n 5

struct coord{
    int i;
    int j;
};
typedef struct coord coord;

double fixdig(double, int);
void triang1(double[][n], int, int); //sem pivotação e pivotação parcial
void triang2(double[][n], int, coord); //com pivotação total
void soluc(double[][n], double[]);
int busca_pivo(double[][n], int);
coord busca_Pivo(double[][n], int);
void troca_linha(double[][n], int, int, int);
double abs_mod(double);
#endif // BIBLIOTECA_SISTEMASLINEARES_H_INCLUDED
