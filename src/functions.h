#ifndef FUNCTION_HEADER
#define FUNCTION_HEADER
#include <armadillo>

using namespace arma;

// diagonalize_armadillo.cpp:
cx_vec diagonalize_arma(mat A);

// create_tridiagonal_matrix.cpp:
void tridag_mat(mat& A, double a, double b, double c, int N);

#endif
