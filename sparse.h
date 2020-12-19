struct sparse {
  int col;
  double val;
};

typedef struct sparse Sparse;

Sparse** sparseCria(int n);
Sparse** criaMatrizA(int n);
void sparseMultmv (int n, Sparse** A, double* v, double* w);
double sparseGet(int i, int j, Sparse** A);
Sparse** sparseMultm(int n, Sparse** A, Sparse** B);