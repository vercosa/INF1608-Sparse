// metodo dos gradientes conjugados
// Ax = b --> (M^-1)Ax = (M^-1)b

// 1. metodo sem pre-condicionador: M = I

// 2.  metodo com pre-condicionador de Gauss-Seidel: M = (D + wL)D^1(D + wU), com w = 1.0

// 3. pre-condicionador SSOR (symmetric successive over-relaxation)
// A = L + D + U
// M = (D + w*L)*(D^-1)*(D+w*U) com w > 1.0.

#include <stdio.h>
#include <math.h>
#include "sparse.h"
#include "matriz.h"
#include "gradconj.h"


double* criaVetorSolucao(int n)
{
    double* x = criavet (n);
    
    for(int i = 0; i < n; i++)
        x[i] = 1.0;
        
    return x;
}

void testa(int n, Sparse** A, double* b, double* xbarra, double* xsol) {
  int i, iter;
  double dif = 0.0;

  iter = GradConj(n, A, b, xbarra, 10e-7);

  for(i = 0; i < n; i++) {
    //   printf("%f\n\n",xbarra[i]);
    dif += fabs((xbarra[i] - xsol[i])/xsol[i]);
  }
  dif = dif / n;
  printf("Erro: %g%%\n", dif*100);
  printf("Iteracoes: %d\n", iter);
  printf("\n");
}

int main ()
{
  int n = 100;
  Sparse** A;
  double* x;
  double* xBarra;
  double* b;
  
  printf("Metodo sem pre-condicionador\n");
  A = criaMatrizA(n);
  x = criaVetorSolucao(n);
  b = criavet(n);
  xBarra = criavet(n);
  
    // for (int i = 0; i < n; i++)
    // {
    //       printf("i = %d --- A[i][0].col = %d A[i][0].val= %.2f\n", i, A[i][0].col, A[i][0].val);
    //     printf("A[i][1].col = %d A[i][1].val= %.2f\n",A[i][1].col, A[i][1].val);
    //     printf("A[i][2].col = %d A[i][2].val= %.2f\n",  A[i][2].col, A[i][2].val);
    //     printf("A[i][3].col = %d A[i][3].val= %.2f\n",  A[i][3].col, A[i][3].val);
    //     printf("A[i][4].col = %d A[i][4].val= %.2f\n\n",  A[i][4].col, A[i][4].val);
    // }
  sparseMultmv(n, A, x, b);
//   printf("%f %f %f %f %f %f \n", b[0],  b[1], b[2], b[3], b[4], b[5]);
      testa(n, A, b, xBarra, x);
  return 0;
}
