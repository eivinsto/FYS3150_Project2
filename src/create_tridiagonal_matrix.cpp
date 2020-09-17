#include <armadillo>
#include "functions.h"

using namespace arma;

void tridag_mat(mat& A, double a, double b, double c, int N) {
  A(0,0) = b;
  A(0,1) = c;
  A(N-1,N-1) = b;
  A(N-1,N-2) = a;
  for (int i = 1; i<N-1; i++){
    A(i,i-1) = a;
    A(i,i) = b;
    A(i,i+1) = c;
  }
}
