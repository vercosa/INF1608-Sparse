// metodo dos gradientes conjugados
// Ax = b --> (M^-1)Ax = (M^-1)b

// 1. metodo sem pre-condicionador: M = I

// 2.  metodo com pre-condicionador de Gauss-Seidel: M = (D + wL)D^1(D + wU), com w = 1.0

// 3. pre-condicionador SSOR (symmetric successive over-relaxation)
// A = L + D + U
// M = (D + w*L)*(D^-1)*(D+w*U) com w > 1.0.

#include <stdio.h>

int main ()
{
  printf ("Hello World");

  return 0;
}
