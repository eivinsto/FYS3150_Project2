#ifndef FUNCTION_HEADER
#define FUNCTION_HEADER
#include <armadillo>

using namespace arma;

// analytic_eigen_val.cpp:
void anal_eig(vec& eigval, double a, double d, int N);

// functions.cpp
void diagonalize_arma(vec& eigval, mat A);
void tridag_mat(mat& A, double a, double d, int N);
void qdot_matrix(mat& A, double rho_max, int N);
void qdot_matrix_double(mat& A, double rho_max, double omega_r, int N);

#endif
