#include <iostream>
#include <armadillo>
#include "jacobi_solver.h"
#include "functions.h"


using namespace arma;
using namespace std;

int main(int argc, char const *argv[])
{
  int N = atoi(argv[1]);
  double rho_max = 10;
  mat A = zeros<mat>(N,N);
  mat R = eye<mat>(N,N);

  qdot_matrix(A,rho_max,N);

  jacobi_solver jacobi(A, R, N);
  jacobi.solve();

  cout << sort(A.diag());
  return 0;
}
