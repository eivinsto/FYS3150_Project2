#ifndef FUNCTION_HEADER
#define FUNCTION_HEADER
#include <armadillo>

using namespace arma;

// diagonalize_armadillo.cpp:
vec diagonalize_arma(mat A);

// create_tridiagonal_matrix.cpp:
void tridag_mat(mat& A, double a, double b, double c, int N);

// analytic_eigen_val.cpp:
vec anal_eig(double d, double a, int N);

#endif
