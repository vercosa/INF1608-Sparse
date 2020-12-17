// int GradConj (int n, double** A, double* b, double* x, double tol);
int GradConj (int n, Sparse** A, double* b, double* x, double tol);

Sparse** precondSsor(int n, Sparse** A, double* b, double w);