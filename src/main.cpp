#include <iostream>
#include <armadillo>
#include <string>
#include "jacobi_solver.h"
#include "functions.h"


using namespace arma;
using namespace std;

int main(int argc, char const *argv[]) {
  int N = atoi(argv[1]);
  string strarg = string(argv[2]);

  if (strarg == "Plot") {
    double pmin = atof(argv[3]), pmax = atof(argv[4]);
    double h = (pmax-pmin)/(double(N));
    double a = -1/(h*h), d = 2/(h*h);
    mat A = zeros<mat>(N,N);
    vec anal_eigvals(N);
    vec comp_eigvals(N);

    //initialize matrices and vector
    tridag_mat(A, a, d, N);
    anal_eig(anal_eigvals, a, d, N);
  }
  return 0;
}

// int main(int argc, char const *argv[])
// {
//   int N = atoi(argv[1]);
//   double rho_max = 10;
//   mat A = zeros<mat>(N,N);
//   mat R = eye<mat>(N,N);
//
//   qdot_matrix_double(A,rho_max,0.01,N);
//
//   jacobi_solver jacobi(A, R, N);
//   jacobi.solve();
//
//   sort(A.diag()).print();
//   return 0;
// }
