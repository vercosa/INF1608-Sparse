#include "sparse.h"
#include <stdlib.h>
#include <stdio.h>


Sparse** sparseCria(int n)
{
  Sparse** A = malloc(sizeof(Sparse*) * n);

  for (int i = 0; i < n; i++)
  {
    A[i] = malloc(sizeof(Sparse));
    A[i][0].col = -1;
    
  }

  return A;
}

Sparse** criaMatrizA(int n)
{
    Sparse** A = sparseCria(n);
    
    for (int i = 0; i < n; i++)
    {
        A[i] = malloc(sizeof(Sparse) * 5);
        
        A[i][0].col = i;
        A[i][0].val = i+1;
        
        // checar se i+1 > n
        if (i+1 < n)
        {
            A[i][1].col = i+1;
            A[i][1].val = 0.5;
        }
        else
        {
            A[i][1].col = -1;
            continue;
        }
        
        // checar se i+2 > n
        if( i+2 <n )
        {
            A[i][2].col = i+2;
            A[i][2].val = 0.5;
        }
        else
        {
            A[i][2].col = -1;
            continue;
        }
        
        if (i != 0)
        {
            if(2*i < n)
            {
                // checar se 2*i é igual a i+1 ou i+2
                if (2*i != i+1 && 2*i != i+2)
                {
                    A[i][3].col = 2*i;
                    A[i][3].val = 0.5;
                }
                else
                {
                    A[i][3].col = -1;
                    continue;                    
                }
            }
            else
            {
                A[i][3].col = -1;
                continue;
            }
        }
        else
        {
            A[i][3].col = -1;
            continue;
        }
        
        A[i][4].col = -1;
        
        // printf("i = %d --- A[i][0].col = %d A[i][0].val= %.2f\n", i, A[i][0].col, A[i][0].val);
        // printf("A[i][1].col = %d A[i][1].val= %.2f\n",A[i][1].col, A[i][1].val);
        // printf("A[i][2].col = %d A[i][2].val= %.2f\n",  A[i][2].col, A[i][2].val);
        // printf("A[i][3].col = %d A[i][3].val= %.2f\n",  A[i][3].col, A[i][3].val);
        // printf("A[i][4].col = %d A[i][4].val= %.2f\n\n",  A[i][4].col, A[i][4].val);
    }
    return A;
}

void sparseMultmv (int n, Sparse** A, double* v, double* w)
{
    int i, j, col;
    
    for (i = 0; i < n; i++)
    {
        w[i] = 0.0;
        
        for (j = 0; j < n; j++)
        {

            col = A[i][j].col;
            
            if (col == -1)
                break;

            if (j == 0 && col > 0)
            {
                // soma os itens abaixo da diagonal
                int k = col-1;
                while (k >= 0)
                {
                    // k vai ser a linha da matriz superior
                    for(int p = 0; p < 5; p++)
                    {
                        if (A[k][p].col == i)
                            w[i] += A[k][p].val * v[col];
                    }
                    k--;
                }
            }

            w[i] += A[i][j].val * v[col];

        }
    }
}
double sparseGet(int i, int j, Sparse** A) {
  int k = 0;

  while(A[i][k].col != -1) {
    if(A[i][k].col == j) {
      return A[i][k].val;
    }
    k++;
  }
  return 0;
}
Sparse** sparseMultm(int n, Sparse** A, Sparse** B) {
  int i, j, k, m;
  double num;
  Sparse** C = sparseCria(n);

  for (i = 0; i < n; i++) {
    m = 0;
    C[i] = malloc(sizeof(Sparse) * (n+1));
    for (k = 0; k < n; k++) {
      num = 0.0;

      for (j = 0; j < n; j++) {
        if(A[i][j].col == -1) break;
        num += A[i][j].val * sparseGet(A[i][j].col, k, B);
      }

      if(num != 0) {
        C[i][m].col = k;
        C[i][m].val = num;
        m++;
      }
    }
    C[i][m].col = -1;
  }

  return C;
}