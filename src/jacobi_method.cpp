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


}

void max_off_diag(mat& A, int& k, int& l, int N)
{
  
}
