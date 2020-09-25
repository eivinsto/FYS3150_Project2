#include <iostream>
#include <armadillo>
#include <string>
#include "jacobi_solver.h"
#include "functions.h"


using namespace arma;
using namespace std;

int main(int argc, char const *argv[]) {
  int N = atoi(argv[1]);
  string strarg = string(argv[2]);

  if (strarg == "toeplitz") {
    double pmin = atof(argv[3]), pmax = atof(argv[4]);
    double h = (pmax-pmin)/(double(N));
    double a = -1/(h*h), d = 2/(h*h);
    mat A = zeros<mat>(N,N);
    mat R(N,N);
    mat anal_eigvec = zeros<mat>(N,N);
    vec anal_eigvals(N);

    //initialize matrices and vector
    tridag_mat(A, a, d, N);
    anal_eig(anal_eigvals, anal_eigvec, a, d, N);
    uvec anal_inx = stable_sort_index(anal_eigvals);

    // jacobi_solver
    jacobi_solver jacobi(A, R, N);
    jacobi.solve();
    uvec comp_indx = stable_sort_index(A.diag());

    cout << (A.diag())[comp_indx[0]] << endl;
    cout << anal_eigvals[anal_inx[0]] << endl;
    (anal_eigvec.col(anal_inx[2])).print();
    cout << endl;
    // mat fac = R.t()/anal_eigvec;
    (R.col(comp_indx[2])).print();
    // cout << endl;
    // (A.diag()).print();
    // (anal_eigvec % fac).print();
    cout << norm_dot(R.col(comp_indx[2]), anal_eigvec.col(2)) << endl;
    // anal_eigvals.print();
    // (R.print();
  }
  return 0;
}

// int main(int argc, char const *argv[])
// {
//   int N = atoi(argv[1]);
//   double rho_max = 10;
//   mat A = zeros<mat>(N,N);
//   mat R = eye<mat>(N,N);
//
//   qdot_matrix_double(A,rho_max,0.01,N);
//
//   jacobi_solver jacobi(A, R, N);
//   jacobi.solve();
//
//   sort(A.diag()).print();
//   return 0;
// }
