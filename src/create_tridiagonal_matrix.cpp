#include <armadillo>
#include "functions.h"

using namespace arma;

void tridag_mat(mat& A, double a, double d, int N) {
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
