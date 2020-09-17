#include <iostream>
#include <armadillo>
#include <cmath>
#include "functions.h"

void jacobi_method(mat& A, mat& R, int N)
{
  /* Function that takes an NxN matrix A, an NxN matrix R and an int N, and
  ** uses the jacobi method to calculate the eigenvectors and eigenvalues of A.
  ** The eigenvectors are stored in R, and the returned matrix A has the
  ** eigenvalues as its diagonal elements.
  */

  // Cleaning eigenvector matrix (making sure it is an identity matrix)
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

  int k,l;   // To store the location of the largest element not on the diagonal
  double epsilon = 1.0e-8;  // Tolerance
  double max_number_iterations = double(N)*double(N)*double(N);
  int iterations = 0;

  double max_off_diag = max_offdiag(A, &k, &l, N);

  /* Perform rotations until the tolerance is satisfied, or the max number of
  ** iterations has been reached.
  */
  while ( fabs(max_off_diag)>epsilon && double(iterations)<max_number_iterations){
    max_off_diag = max_offdiag(A, &k, &l, N);
    rotate(A, R, k, l, N);
    iterations++;
  }

}

double max_offdiag(mat& A, int& k, int& l, int N)
{
  /* Function that returns the largest element in A that is not on the diagonal
  ** and stores its indices in k and l.
  */

  double max = 0.0;
  for(int i = 0; i < n; i++ ){
    for(int j = i + 1; j < n; j++ ){
      if( fabs(A[i][j]) > max ){
        max = fabs(A[i][j]);
        *l = i;
        *k = j;
      }
    }
  }
  return max;
}

void rotate(mat& A, mat& R, int k, int l, int N)
{
  /* Function that performs the rotation of A. Modifies A and R.
  */

  double s,c;

  // Find elements of rotation matrix
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
    s = c*t;
  }
  else {
    c = 1.0;
    s = 0.0;
  }

  // Create some variables needed later when overwriting the values stored in A.
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
