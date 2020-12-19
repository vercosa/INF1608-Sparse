struct sparse {
  int col;
  double val;
};

typedef struct sparse Sparse;

Sparse** sparse_cria(int n);
Sparse** criaMatrizA(int n);
void sparseMultmv (int n, Sparse** A, double* v, double* w);
double sparseGet(int i, int j, Sparse** A);