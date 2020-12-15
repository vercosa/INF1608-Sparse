// metodo dos gradientes conjugados
// Ax = b --> (M^-1)Ax = (M^-1)b

// 1. metodo sem pre-condicionador: M = I

// 2.  metodo com pre-condicionador de Gauss-Seidel: M = (D + wL)D^1(D + wU), com w = 1.0

// 3. pre-condicionador SSOR (symmetric successive over-relaxation)
// A = L + D + U
// M = (D + w*L)*(D^-1)*(D+w*U) com w > 1.0.

#include <stdio.h>
#include "sparse.h"
#include "matriz.h"

double* criaVetorSolucao(int n)
{
    double* x = criavet (n);
    
    for(int i = 0; i < n; i++)
        x[i] = 1.0;
        
    return x;
}

int main ()
{
  int n = 6;
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
  printf("\n\n%f %f %f %f %f %f ", b[0],  b[1], b[2], b[3], b[4], b[5]);

  return 0;
}
