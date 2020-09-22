#include "catch.hpp"
#include "jacobi_solver.h"
#include "functions.h"


TEST_CASE("Testing eigenvalues with armadillo diagonalizer."){
    int n = 3;
    double pmin = 0, pmax = 1;
    double h = (pmax-pmin)/(double(n));
    double a = -1/(h*h), d = 2/(h*h);
    mat A = zeros<mat>(n,n);
    vec anal_eigvals(n);
    vec comp_eigvals(n);

    //initialize matrices and vector
    tridag_mat(A, a, d, n);
    anal_eig(anal_eigvals, a, d, n);

    // finding eigenvalues with armadillo
    eig_sym(comp_eigvals, A);

    anal_eigvals = sort(anal_eigvals);
    comp_eigvals = sort(comp_eigvals);

    // testing results
    for (int i = 0; i < n; i++) {
      REQUIRE(comp_eigvals[i] == Approx(anal_eigvals[i]).epsilon(0.001));
    }
}


TEST_CASE("Testing eigenvalues with jacobi_solver."){
    int n = 3;
    double pmin = 0, pmax = 1;
    double h = (pmax-pmin)/(double(n));
    double a = -1/(h*h), d = 2/(h*h);
    mat A = zeros<mat>(n,n);
    mat R = eye<mat>(n,n);
    vec anal_eigvals(n);

    //initialize matrices and vector
    tridag_mat(A, a, d, n);
    anal_eig(anal_eigvals, a, d, n);

    // finding eigenvalues with armadillo
    jacobi_solver jacobi(A, R, n);
    jacobi.solve();

    anal_eigvals = sort(anal_eigvals);
    vec comp_eigvals = sort(A.diag());

    // testing results
    for (int i = 0; i < n; i++) {
      REQUIRE(comp_eigvals[i] == Approx(anal_eigvals[i]).epsilon(0.001));
    }
}


TEST_CASE("Testing orthogonality of eigenvectors from jacobi_solver."){
    int n = 3;
    double tol = 1e-14; // tolerance for comparison to zero
    double pmax = 1;
    double h = pmax/(double(n));
    double a = -1/(h*h), d = 2/(h*h);

    //initialize matrix
    mat A = zeros<mat>(n,n);
    mat R = eye<mat>(n,n);
    tridag_mat(A, a, d, n);

    // finding eigenvectors with jacobi_solver
    jacobi_solver jacobi(A, R, n);
    jacobi.solve();

    // taking dot-product between orthogonal eigenvectors
    double result1 = abs(dot( R.row(0), R.row(1) ));
    double result2 = abs(dot( R.row(0), R.row(2) ));
    double result3 = abs(dot( R.row(1), R.row(2) ));

    // checking if dot-product of eigenvectors is sufficently close to zero
    REQUIRE(result1 < tol);
    REQUIRE(result2 < tol);
    REQUIRE(result3 < tol);
}


TEST_CASE("Testing  jacobi_solver."){
    int k, l;
    int n = 5;
    int p = 3, q = 4;
    // double pmax = 1;
    mat A = eye<mat>(n,n);
    mat R = eye<mat>(n,n);

    // setting largest non-zero diagonal at known indices
    A(p, q) = -8.0;

    // finding eigenvalues with armadillo
    jacobi_functions jacobi(A, R, n);
    jacobi.max_offdiag(&k, &l);
    REQUIRE(l == p);
    REQUIRE(k == q);
    REQUIRE(A(l, k) == Approx(-8.0).epsilon(0.001));

}
