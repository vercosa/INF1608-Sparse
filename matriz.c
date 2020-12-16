#include "matriz.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* criavet (int n)
{
    if ( n<0 )
        return NULL;
    
    return (double*)malloc(n * sizeof(double));
}

void liberavet (double* v)
{
    free(v);
}

double prodescalar (int n, const double* v, const double* w)
{
    double result = 0.0;
    int i;
    
    if ( v==NULL  || w==NULL )
        return result;
        
    for (i = 0; i<n; i++)
        result = result + (v[i] * w[i]);
        
    return result;
}

double norma2 (int n,const double* v)
{
    double result = 0.0;
    int i;
    
    for (i=0; i<n; i++)
        result = result + pow(v[i], 2);
        
    return sqrt(result);
}

void multvs (int n, double* v, double s, double *w)
{
    int i;
    
    if ( v==NULL  || w==NULL )
        return;
        
    for (i = 0; i<n; i++)
        w[i] = v[i] * s;
}

double** criamat (int m, int n)
{
    double** mat;
    int i;
    
    if (m < 0)
        return NULL;
        
    mat = (double**)malloc(m * sizeof(double*));
    
    if (mat == NULL)
        return mat;
        
    for(i=0; i < m; i++)
        mat[i] = criavet (n);
        
    return mat;
}

double** criamattri (int m)
{
    double** mat;
    int i;
    
    if (m < 0)
        return NULL;
        
    mat = (double**)malloc(m * sizeof(double*));
    
    if (mat == NULL)
        return mat;
        
    for(i=0; i < m; i++)
        mat[i] = criavet (i+1);
        
    return mat;
}

void liberamat (int m, double** A)
{
    int i;
    
    for(i = 0; i <m ; i++)
        liberavet(A[i]);
        
    free(A);
}

void transposta (int m, int n, double** A, double** T)
{
    int i, j;
    
    for(i=0; i < m; i++)
        for(j = 0; j < n; j++)
            T[j][i] = A[i][j];
}

void multmv (int m, int n, double** A, double* v, double* w)
{
    int i;
    
    if (A == NULL  || v == NULL  || w == NULL)
        return;
        
    for(i=0; i < m; i++)
        w[i] = prodescalar(n, A[i], v);
}

void multmm (int m, int n, int q, double** A, double** B, double** C)
{
    double** B_transposta =criamat(q, n);
    transposta(n, q, B, B_transposta);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < q; j++) {
            C[i][j] = prodescalar(n, A[i], B_transposta[j]);
        }
    }
}

double vet_dot(int n, double* v, double* w)
{
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += v[i] * w[i];
  }
  return s;
}

void vet_mults(int n, double* v, double s, double* w)
{
  int i;
  for (i = 0; i < n; i++) {
    w[i] = s * v[i];
  }
}

void vet_soma(int n, double* a, double* b, double* c)
{
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] + b[i];
  }
}

void vet_subtrai(int n, double* a, double* b, double* c)
{
  int i;
  for (i = 0; i < n; i++) {
    c[i] = a[i] - b[i];
  }
}

void vet_copia(int n, double* src, double* dst) {
  int i;
  for (i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}