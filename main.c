// metodo dos gradientes conjugados
// Ax = b --> (M^-1)Ax = (M^-1)b

// 1. metodo sem pre-condicionador: M = I

// 2.  metodo com pre-condicionador de Gauss-Seidel: M = (D + wL)D^1(D + wU), com w = 1.0

// 3. pre-condicionador SSOR (symmetric successive over-relaxation)
// A = L + D + U
// M = (D + w*L)*(D^-1)*(D+w*U) com w > 1.0.

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "sparse.h"
#include "matriz.h"
#include "gradconj.h"

double *criaVetorSolucao(int n)
{
  double *x = criavet(n);

  for (int i = 0; i < n; i++)
    x[i] = 1.0;

  return x;
}

void testa(int n, Sparse **A, double *b, double *xbarra, double *xsol, int testType, int w)
{
  int i, iter;
  double dif = 0.0;
  clock_t end, start;

  start = clock();

  if (testType == 2){
    A = precondGauss(n, A, b);
  }
  if (testType == 3){
    A = precondSsor(n, A, b, w);
  }

  iter = GradConj(n, A, b, xbarra, 10e-7);

  for (i = 0; i < n; i++)
    dif += fabs((xbarra[i] - xsol[i]) / xsol[i]);

  dif = dif / n;

  end = clock();
  printf("Erro: %g%%\n", dif * 100);
  printf("Iteracoes: %d\n", iter);
  printf("Durou %lf segundos\n\n", (double)(end - start) / CLOCKS_PER_SEC);
  printf("\n");
}

int main()
{
  int n = 10000;
  Sparse **A;
  double *x;
  double *xBarra;
  double *b;
  float w = 1.1;

  printf("Metodo sem pre-condicionador\n");
  A = criaMatrizA(n);
  x = criaVetorSolucao(n);
  b = criavet(n);
  xBarra = criavet(n);
  sparseMultmv(n, A, x, b);
  testa(n, A, b, xBarra, x, 1, 0);

  printf("Metodo de Gauss\n");
  A = criaMatrizA(n);
  x = criaVetorSolucao(n);
  b = criavet(n);
  xBarra = criavet(n);
  sparseMultmv(n, A, x, b);
  testa(n, A, b, xBarra, x, 2, 0);

  while (w < 2.1)
  {
    printf("Metodo de SSOR, w = %f\n", w);
    A = criaMatrizA(n);
    x = criaVetorSolucao(n);
    b = criavet(n);
    xBarra = criavet(n);
    sparseMultmv(n, A, x, b);
    testa(n, A, b, xBarra, x, 3, w);
    w += 0.1;
  }

  printf("-----------------------------\n\n");

  return 0;
}
