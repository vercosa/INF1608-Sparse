#include "matriz.h"
#include "sparse.h"
#include <stdio.h>

int GradConj (int n, Sparse** A, double* b, double* x, double tol)
{
    int count = 0;
    // double rTr = 0.0;
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
    // multmv(n,n,A,x,Ax);
    
    vet_copia(n, r, d);
    
    // for(int i=0; i<n; i++)
    //     rTr = rTr + (r[i]*r[i]);

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
        printf("%d\n", count);
    }

    return count;
}

// int GradConjPreCond (int n, double** A, double* b, double* x, double tol)
// {
//     int count = 0;
//     double rTz = 0.0;
//     double alfa, beta;
//     double *d = criavet(n);
//     double *r = criavet(n);
//     double *z = criavet(n);
//     double* Ax = criavet(n);
    
//     multmv(n,n,A,x,Ax);
    
//     for(int i=0; i<n; i++)
//     {
//         r[i]=b[i]-Ax[i];
//         d[i]=r[i]/A[i][i];
//         z[i]=d[i];
//     }
    
//     for(int i=0; i<n; i++)
//         rTz = rTz + (r[i]*z[i]);

//     for(int i=0; i<n; i++)
//     {
//         double dTAd = 0.0;
//         double rTzProx = 0.0;
        
//         if(norma2(n,r)<tol)
//         {
//             break;
//         }
    
//         multmv(n,n,A,d,Ax);
        
//         for(int k=0; k<n; k++)
//             dTAd = dTAd + (d[k]*Ax[k]);

//         alfa = rTz/dTAd;

//         for(int k=0; k<n; k++)
//         {
//             x[k] = x[k] + alfa * d[k];
//             r[k] = r[k] - alfa * Ax[k];
//             z[k] = z[k]/A[k][k];
//         }

//         for(int k=0; k<n; k++)
//         {
//             rTzProx = rTzProx + r[k]*z[k];
//         }

//         beta = rTzProx / rTz;

//         for(int k=0; k<n; k++)
//             d[k] = z[k] + beta * d[k];

//         rTz = rTzProx;
//         count++;
//     }

//     liberavet(r);
//     liberavet(d);
//     liberavet(z);
//     liberavet(Ax);
    
//     return count;
// }