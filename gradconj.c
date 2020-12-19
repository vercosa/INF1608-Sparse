#include "matriz.h"
#include "sparse.h"
#include <stdio.h>
#include <stdlib.h>

int GradConj (int n, Sparse** A, double* b, double* x, double tol)
{
    int count = 0;

    double alfa, beta;
    double *d = criavet(n);
    double *r = criavet(n);
    
    double* Ad = criavet(n);
    double* ax = criavet(n);
    
    double* xn = criavet(n);
    double* ad = criavet(n);
    
    double* rn = criavet(n);
    double* aAd = criavet(n);
    
    double* dn = criavet(n);
    double* bd = criavet(n);
    
    sparseMultmv(n, A, x, ax);
    vet_subtrai(n, b, ax, r);
    
    vet_copia(n, r, d);
    
    for(int i=0; i<n; i++)
    {
        if(norma2(n,r)<tol)
            break;
            
        // alfa_k = rTk . rK / dTk . A dk
        sparseMultmv(n, A, d, Ad);
        alfa = vet_dot(n, r, r) / vet_dot(n, d, Ad);

        // x_{k+1} = xk + alfa_k dk;
        vet_mults(n, d, alfa, ad);
        vet_soma(n, x, ad, xn);
        
        // r_{k+1} = rk - alfa_k A dk;
        vet_mults(n, Ad, alfa, aAd);
        vet_subtrai(n, r, aAd, rn);
        
        // Beta_k = rTn rn / 
        beta = vet_dot(n, rn, rn) / vet_dot(n, r, r);
    
        // d_{k+1} = r_{k+1} + Beta_k dk;
        vet_mults(n, d, beta, bd);
        vet_soma(n, rn, bd, dn);
    
          
        vet_copia(n, xn, x);
        vet_copia(n, rn, r);
        vet_copia(n, dn, d);
        count++;
    }

    return count;
}

Sparse** precondSsor(int n, Sparse** A, double* b, double w)
{
     int i, j, k;
    
    double* b_ = criavet(n);
    
    Sparse** C;
    Sparse** M1 = sparseCria(n);
    Sparse** M2 = sparseCria(n);
    Sparse** M3 = sparseCria(n);
    Sparse** Ma;
    Sparse** Mb;
    Sparse** Minv;
    
    for(i = 0; i < n; i++) {
        M1[i] = malloc(sizeof(Sparse) * n);
        
        M1[i][0].col = i;
        M1[i][0].val = sparseGet(i, i, A);
        
        k = 1;
        for(j = 0; j < n; j++) {
            if(A[i][j].col == -1) break;
            if(A[i][j].col < i) {
                M1[i][k].col = A[i][j].col;
                M1[i][k].val = A[i][j].val * w;
                k++;
            }
        }
        
        M1[i][k].col = -1;
    }
    
    for(i = 0; i < n; i++) {
        M2[i] = malloc(sizeof(Sparse) * 2);
        
        M2[i][0].col = i;
        M2[i][0].val = 1.0/sparseGet(i, i, A);
        
        M2[i][1].col = -1;
    }
    
    for(i = 0; i < n; i++) {
        M3[i] = malloc(sizeof(Sparse) * n);
        
        M3[i][0].col = i;
        M3[i][0].val = sparseGet(i, i, A);
        
        k = 1;
        for(j = 0; j < n; j++) {
            if(A[i][j].col == -1) break;
            if(A[i][j].col > i) {
                M3[i][k].col = A[i][j].col;
                M3[i][k].val = A[i][j].val * w;
                k++;
            }
        }
        
        M3[i][k].col = -1;
    }
    
    Ma = sparseMultm(n, M1, M2);
    Mb = sparseMultm(n, Ma, M3);
    Minv = Mb;
    
    C = sparseMultm(n, Minv, A);
    
    sparseMultmv(n, Minv, b, b_);
    vet_copia(n, b_, b);
    
    return C;
}
Sparse** precond_Gauss(int n, Sparse** A, double* b)
{
     int i, j, k;
    
    double* b_ = criavet(n);
    
    Sparse** C;
    Sparse** M1 = sparseCria(n);
    Sparse** M2 = sparseCria(n);
    Sparse** M3 = sparseCria(n);
    Sparse** Ma;
    Sparse** Mb;
    Sparse** Minv;
    
    for(i = 0; i < n; i++) {
        M1[i] = malloc(sizeof(Sparse) * n);
        
        M1[i][0].col = i;
        M1[i][0].val = sparseGet(i, i, A);
        
        k = 1;
        for(j = 0; j < n; j++) {
            if(A[i][j].col == -1) break;
            if(A[i][j].col < i) {
                M1[i][k].col = A[i][j].col;
                M1[i][k].val = A[i][j].val;
                k++;
            }
        }
        
        M1[i][k].col = -1;
    }
    
    for(i = 0; i < n; i++) {
        M2[i] = malloc(sizeof(Sparse) * 2);
        
        M2[i][0].col = i;
        M2[i][0].val = 1.0/sparseGet(i, i, A);
        
        M2[i][1].col = -1;
    }
    
    for(i = 0; i < n; i++) {
        M3[i] = malloc(sizeof(Sparse) * n);
        
        M3[i][0].col = i;
        M3[i][0].val = sparseGet(i, i, A);
        
        k = 1;
        for(j = 0; j < n; j++) {
            if(A[i][j].col == -1) break;
            if(A[i][j].col > i) {
                M3[i][k].col = A[i][j].col;
                M3[i][k].val = A[i][j].val ;
                k++;
            }
        }
        
        M3[i][k].col = -1;
    }
    
    Ma = sparseMultm(n, M1, M2);
    Mb = sparseMultm(n, Ma, M3);
    Minv = Mb;
    
    C = sparseMultm(n, Minv, A);
    
    sparseMultmv(n, Minv, b, b_);
    vet_copia(n, b_, b);
    
    return C;
}
