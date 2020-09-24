#ifndef JACOBI_SOLVER_H
#define JACOBI_SOLVER_H

#include <armadillo>

using namespace arma;

/* Splitting this into two classes allows us to test the individual functions
** while still keeping them private in terms of the jacobi_solver class.
*/

class jacobi_functions{
/* Class containing some of the functions used in jacobi_solver.
*/
public:
  mat* A;
  mat* R;
  int N,k,l;

  jacobi_functions(mat& Am, mat& Rm, int Nm){
    A = &Am;
    R = &Rm;
    N = Nm;
  }

  double max_offdiag();
  void rotate();

};

class jacobi_solver: private jacobi_functions{
/* Class that finds eigenvalues and eigenvectors of a symmetric NxN matix A.
** After running solve the matrix A will have the eigenvalues along its diagonal
** and R will contain the eigenvectors as rows.
*/
public:
  jacobi_solver(mat& Am, mat& Rm, int Nm)
  : jacobi_functions(Am, Rm, Nm)
  {
  }

  void solve();
};

#endif
