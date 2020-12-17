#include "matriz.h"
#include "sparse.h"
#include <stdio.h>

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

Sparse** precond_ssor(int n, Sparse** A, double* b, double w)
{
    
}
