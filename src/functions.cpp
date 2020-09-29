#include <armadillo>
#include <cmath>
#include "functions.h"

using namespace arma;

void tridag_mat(mat& A, double a, double d, int N) {
  /*
  ** Populates tridiagonal matrix.
  ** Args:
  **   A - NxN - zero-matrix to populate.
  **   a - value of non-zero, off-diagonal elements.
  **   d - value of diagonal elements.
  **   N - dimensionality of matrix.
  */
  A = zeros<mat>(N,N);

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
  /*
  ** calculates analytic eigenvalues for Toeplitz-matrix.
  ** Args:
  **   eigval - vector to populate.
  **   a - value of non-zero, off-diagonal elements.
  **   d - value of diagonal elements.
  **   N - dimensionality of matrix.
  */

  // creating vector of j-values
  vec j_val = zeros<vec>(N);
  for (int i = 0; i <= N; i++) {
    j_val[i] = 1 + i;
  }

  // calculating eigenvalues element-wise
  eigval = d + 2*a*cos(j_val*M_PI/(N+1));
}

void anal_eig(vec& eigval, mat& eigvec, double a, double d, int N) {
  /*
  ** calculates analytic eigenvalues and -vectors for Toeplitz-matrix.
  ** Args:
  **   eigval - vector to populate with eigenvalues.
  **   eigvec - matrix to store eigenvectors in columns of.
  **   a - value of non-zero, off-diagonal elements.
  **   d - value of diagonal elements.
  **   N - dimensionality of matrix.
  */

  // creating vector of j-values
  vec j_val = zeros<vec>(N);
  for (int i = 0; i < N; i++) {
    j_val[i] = 1 + i;
  }

  // calculating eigenvalues element-wise
  eigval = d + 2*a*cos(j_val*M_PI/(N+1));

  // storing eigenvectors as columns in eigvec matrix
  for (int i = 0; i < N; i++) {
    eigvec.col(i) = sin((i+1)*j_val*M_PI/(N+1));
  }
}

void qdot_matrix(mat& A, double rho_max, int N){
  /*
  ** populate tridag matrix for single-electron problem.
  ** Args:
  **   A - matix to populate.
  **   rho_max - maximum value of rho.
  **   N - dimensionality of matrix.
  */

  A = zeros<mat>(N,N);
  double h = rho_max/(N+1);
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

void qdot_matrix_double(mat& A, double rho_max, double omega_r, int N){
  /*
  ** populate tridag matrix for two-electron problem.
  ** Args:
  **   A - matix to populate.
  **   rho_max - maximum value of rho.
  **   omega_r - value of omega_r.
  **   N - dimensionality of matrix.
  */

  A = zeros<mat>(N,N);
  double h = rho_max/(N+1);
  double hh = h*h;
  double e = -1/hh;
  double d0 = 2/hh;
  double wr2 = omega_r*omega_r;
  vec rho = linspace(h,rho_max-h,N);

  A(0,0) = d0 + rho(0)*rho(0);
  for (int i = 1; i<N; ++i){
    A(i-1,i) = e;
    A(i,i-1) = e;
    A(i,i) = d0 + wr2*rho(i)*rho(i) + 1/rho(i);
  }
}

void qdot_matrix_double_eigval(vec& e, double omega_r, int N){
  /*
  ** calculates analytic approximate of eigenvalues for two-electron problem.
  ** Args:
  **   e - vector to populate.
  **   omega_r - value of omega_r.
  **   N - dimensionality of matrix.
  */
  for (int i = 0; i<N; ++i){
    e(i) = 3*cbrt((omega_r/2)*(omega_r/2)) + sqrt(3)*omega_r*(2*i+1);
  }
}
