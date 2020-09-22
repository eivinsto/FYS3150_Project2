#include <armadillo>
#include <cmath>
#include "functions.h"

using namespace arma;

void tridag_mat(mat& A, double a, double d, int N) {
  // Generate general tridiagonal matrix
  A(0,0) = d;
  A(0,1) = a;
  A(N-1,N-1) = d;
  A(N-1,N-2) = a;
  for (int i = 1; i<N-1; i++){
    A(i,i-1) = a;
    A(i,i) = d;
    A(i,i+1) = a;
  }
}

void anal_eig(vec& eigval, double a, double d, int N) {
  // Analytic eigen values of general tridiagonal matrix

  // creating vector of j-values
  vec j_val = zeros<vec>(N);
  for (int i = 0; i < N + 1; i++) {
    j_val[i] = 1 + i;
  }

  // calculating eigenvalues element-wise
  eigval = d - 2*a*cos(j_val*M_PI/(N+1));
}

void qdot_matrix(mat& A, double rho_max, int N){
  // Generates matrix for quantum dot problem
  double h = rho_max/double(N);
  double hh = h*h;
  double e = -1/hh;
  double d0 = 2/hh;
  vec rho = linspace(h,rho_max-h,N);

  A(0,0) = d0 + rho(0)*rho(0);
  for (int i = 1; i<N; ++i){
    A(i-1,i) = e;
    A(i,i-1) = e;
    A(i,i) = d0 + rho(i)*rho(i);
  }
}
