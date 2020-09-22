#ifndef FUNCTION_HEADER
#define FUNCTION_HEADER
#include <armadillo>

using namespace arma;

// diagonalize_armadillo.cpp:
void diagonalize_arma(vec& eigval, mat A);

// create_tridiagonal_matrix.cpp:
void tridag_mat(mat& A, double a, double d, int N);

// analytic_eigen_val.cpp:
void anal_eig(vec& eigval, double d, double a, int N);

#endif
