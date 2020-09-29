#include <iostream>
#include <armadillo>
#include <string>
#include <time.h>
#include "jacobi_solver.h"
#include "functions.h"


using namespace arma;
using namespace std;

int main(int argc, char const *argv[]) {
  /*
  ** This program solves the problems from the report.
  ** For best user experience avoid using this program directly.
  ** Instead, use the python script in the root folder of the repo.
  */
  int N = atoi(argv[1]) - 1;  // retriveing N-1 as N.
  string strarg = string(argv[2]);  // retriveing which problem to solve.


  if (strarg == "toeplitz" || strarg == "to") {
    /*
    ** Solving buckling beam problem from report.
    ** Results are written to files.
    */

    // retriveing rho_min and rho_max:
    double pmin = atof(argv[3]), pmax = atof(argv[4]);
    double h = (pmax-pmin)/(N+1);  // calculating steplength
    double a = -1/(h*h), d = 2/(h*h);  // calculating non-zero matix elements
    mat A(N,N);  // declaring toeplitz matrix
    mat R(N,N);  // matrix for storing eigenvectors
    vec anal_eigvals(N);  // vector for analytic eigenvalues
    mat anal_eigvecs(N,N);  // matrix for analytic eigenvectors

    //initialize matrices and vector
    tridag_mat(A, a, d, N);  // populating toeplitz matrix
    anal_eig(anal_eigvals, anal_eigvecs, a, d, N);  // calculating analytic solution

    // saving analytic solution
    anal_eigvals.save("anal_eigvals_" + to_string(N+1) + ".dat", arma_ascii);
    anal_eigvecs.save("anal_eigvecs_" + to_string(N+1) + ".dat", arma_ascii);

    // solving with jacobi
    jacobi_solver jacobi(A, R, N);
    jacobi.solve();

    // saving numerical result to files
    A.save("comp_eigvals_" + to_string(N+1) + ".dat", arma_ascii);
    R.save("comp_eigvecs_" + to_string(N+1) + ".dat", arma_ascii);
  }


  if (strarg == "single" || strarg == "s") {
    /*
    ** Solving single-electron problem from report.
    ** Results are written to file.
    */
    double pmax = atof(argv[3]);
    mat A(N,N);
    mat R(N,N);

    // initialize matrices and vector
    qdot_matrix(A, pmax, N);

    // calculates numerical result
    jacobi_solver jacobi(A, R, N);
    jacobi.solve();

    vec eigvals = A.diag();  // extracts numerical eigenvalues from A.

    // saves eigenvalues to file
    eigvals.save("quantum_eigvals_" + to_string(N+1) + ".dat", arma_ascii);
  }


  if (strarg == "double" || strarg == "d") {
    /*
    ** Solving two-electron problem from report.
    ** Results are written to files.
    */
    double pmax = atof(argv[3]);  // retriveing rho_max
    double omega_r = atof(argv[4]);  // retriveing omega_r
    mat A(N,N);
    mat R(N,N);
    vec anal_eigvals(N);

    // calculates analytic approximate energies
    qdot_matrix_double_eigval(anal_eigvals, omega_r, N);

    // initializing tridag matrix A.
    qdot_matrix_double(A, pmax, omega_r, N);

    // solving numerical energies
    jacobi_solver jacobi(A, R, N);
    jacobi.solve();

    vec eigvals = A.diag();  // extracting numerical eigenvalues

    // saving analytic approximate and numerical eigenvalues to files
    anal_eigvals.save("double_qdot_anal_eigvals_" + to_string(N+1) + ".dat", arma_ascii);
    eigvals.save("double_qdot_num_eigvals_" + to_string(N+1) + ".dat", arma_ascii);
  }


  if (strarg == "benchmark" || strarg == "b") {
    /*
    ** Performs benchmark of jacobi solver and Armadillo eig_sym.
    ** Results are written to file.
    */

    // calculating steplength and non-zero matrix elements
    double h = 1/(N+1);
    double a = -1/(h*h), d = 2/(h*h);

    // declaring matrices and vectors for results
    mat A(N,N);
    mat R(N,N);
    mat eigvecs(N,N);
    vec eigvals(N);
    mat benchmark_times(N,2);

    // declaring clock_t objects for timing
    clock_t start, finish;

    // performs solve N times for each solver
    for (int i = 0; i < N; i++) {
      cout << "\tPerforming benchmark: " << 100*(i+1)/N << "%\r";
      cout.flush();

      // initializing matrix
      tridag_mat(A, a, d, N);

      // benchmarking eig_sym
      start = clock();
      eig_sym(eigvals, eigvecs, A);
      finish = clock();

      // saving result to matrix
      benchmark_times(i,0) = double(finish - start)/CLOCKS_PER_SEC;

      // creating instance of jacobi_solver, this cleans R
      jacobi_solver jacobi(A, R, N);

      // benchmarking jacobi_solver
      start = clock();
      jacobi.solve();
      finish = clock();

      // saving result to matrix
      benchmark_times(i,1) = double(finish - start)/CLOCKS_PER_SEC;
    }
    cout << "\n\tBenchmark complete!" << endl;

    // saving results to file
    benchmark_times.save("benchmark_times.dat", arma_ascii);
  }
  return 0;
}
