#include "matriz.h"

int GradConj (int n, double** A, double* b, double* x, double tol)
{
    int count = 0;
    double rTr = 0.0;
    double alfa, beta;
    double *d = criavet(n);
    double *r = criavet(n);
    double* Ax = criavet(n);
    
    multmv(n,n,A,x,Ax);
    
    for(int i=0; i<n; i++)
    {
        d[i] = b[i]-Ax[i];
        r[i] = b[i]-Ax[i];
    }
    
    for(int i=0; i<n; i++)
        rTr = rTr + (r[i]*r[i]);

    for(int i=0; i<n; i++)
    {
        double dTAd = 0.0;
        double rTrProx = 0.0;
        
        if(norma2(n,r)<tol)
            break;

        multmv(n,n,A,d,Ax);
        
        for(int k=0; k<n; k++)
            dTAd = dTAd + (d[k]*Ax[k]);

        alfa = rTr/dTAd;

        for(int k=0; k<n; k++)
        {
            x[k] = x[k] + alfa * d[k];
            r[k] = r[k] - alfa * Ax[k];
        }

        for(int k=0; k<n; k++)
            rTrProx = rTrProx + (r[k]*r[k]);

        beta = rTrProx / rTr;

        for(int k=0; k<n; k++)
            d[k] = r[k] + beta * d[k];

        rTr = rTrProx;
        count++;
    }


    liberavet(d);
    liberavet(r);
    liberavet(Ax);
    return count;
}

int GradConjPreCond (int n, double** A, double* b, double* x, double tol)
{
    int count = 0;
    double rTz = 0.0;
    double alfa, beta;
    double *d = criavet(n);
    double *r = criavet(n);
    double *z = criavet(n);
    double* Ax = criavet(n);
    
    multmv(n,n,A,x,Ax);
    
    for(int i=0; i<n; i++)
    {
        r[i]=b[i]-Ax[i];
        d[i]=r[i]/A[i][i];
        z[i]=d[i];
    }
    
    for(int i=0; i<n; i++)
        rTz = rTz + (r[i]*z[i]);

    for(int i=0; i<n; i++)
    {
        double dTAd = 0.0;
        double rTzProx = 0.0;
        
        if(norma2(n,r)<tol)
        {
            break;
        }
    
        multmv(n,n,A,d,Ax);
        
        for(int k=0; k<n; k++)
            dTAd = dTAd + (d[k]*Ax[k]);

        alfa = rTz/dTAd;

        for(int k=0; k<n; k++)
        {
            x[k] = x[k] + alfa * d[k];
            r[k] = r[k] - alfa * Ax[k];
            z[k] = z[k]/A[k][k];
        }

        for(int k=0; k<n; k++)
        {
            rTzProx = rTzProx + r[k]*z[k];
        }

        beta = rTzProx / rTz;

        for(int k=0; k<n; k++)
            d[k] = z[k] + beta * d[k];

        rTz = rTzProx;
        count++;
    }

    liberavet(r);
    liberavet(d);
    liberavet(z);
    liberavet(Ax);
    
    return count;
}