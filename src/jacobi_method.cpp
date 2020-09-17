#include <iostream>
#include <armadillo>
#include <cmath>
#include "functions.h"

void jacobi_method(mat& A, mat& R, int N)
{
  // Cleaning eigenvector matrix
  for (int i; i<N; ++i){
    for (int j; i<N; ++i){
      if (i==j){
        R[i][j] = 1;
      }
      else{
        R[i][j] = 0;
      }
    }
  }

  int k,l;
  double epsilon = 1.0e-8;
  double max_number_iterations = double(N)*double(N)*double(N);
  int iterations = 0;
  double max_off_diag = max_offdiag(A, &k, &l, N);

  while ( fabs(max_off_diag)>epsilon && double(iterations)<max_number_iterations){
    max_off_diag = max_offdiag(A, &k, &l, N);
    rotate(A, R, k, l, N);
    iterations++;
  }

}

double max_offdiag(mat& A, int& k, int& l, int N)
{
  double max = 0.0;
  for (int i; i<N; ++i){
    for (int j; i<N; ++j){
      if (A[i][j]>max){
        max = A[i][j]
        *k = i
        *l = j
      }
    }
  }
  return max;
}

void rotate(mat& A, mat& R, int k, int l, int N)
{
  double s,c;

  if (A[k][l] != 0.0){
    double t,tau;
    tau = (A[l][l]- A[k][k])/(2*A[k][l]);
    if (tau > 0){
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    }
    else {
      t = -1.0/(-tau + sqrt(1.0 + tau*tau));
    }

    c = 1.0/sqrt(1+t*t);
    s = c*t
  }
  else {
    c = 1.0;
    s = 0.0;
  }

  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A[k][k];
  a_ll = A[l][l];

  // Changing the matrix elements with indices k and l
  A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
  A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
  A[k][l] = 0.0; // hard-coding of the zeros
  A[l][k] = 0.0;

  // Changing the remaining elements
  for (int i; i<N; ++i){
    if (i != k && i != l) {
      a_ik = A[i][k];
      a_il = A[i][l];
      A[i][k] = c*a_ik - s*a_il;
      A[k][i] = A[i][k];
      A[i][l] = c*a_il + s*a_ik;
      A[l][i] = A[i][l];
    }

    // Compute the new eigenvectors
    r_ik = R[i][k];
    r_il = R[i][l];
    R[i][k] = c*r_ik - s*r_il;
    R[i][l] = c*r_il + s*r_ik;
  }
}
