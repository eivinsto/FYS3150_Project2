#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include <armadillo>

using namespace arma;

class jacobi_functions{
public:
  mat* A;
  mat* R;
  int N;

  jacobi_functions(mat& Am, mat& Rm, int Nm){
    A = &Am;
    R = &Rm;
    N = Nm;
  }

  double max_offdiag(int* k, int* l);
  void rotate(int k, int l);

};

class jacobi_method: private jacobi_functions{
public:
  jacobi_method(mat& Am, mat& Rm, int Nm)
  : jacobi_functions(Am, Rm, Nm)
  {
  }

  void solve();
};

#endif
