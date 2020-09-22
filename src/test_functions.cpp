#include "catch.hpp"
#include "jacobi_solver.h"
#include "functions.h"


TEST_CASE("Testing eigenvalues with armadillo diagonalizer."){
    int n=3;
    double pmin=0, pmax=10;
    double h = (pmax-pmin)/(double(n));
    mat A = zeros<mat>(n,n);
    vec anal_eigvals = zeros<vec>(n);
    vec comp_eigvals(n);

    //initialize matrices and vector
    tridag_mat(A, -1.0/(h*h), 2.0/(h*h), n);
    anal_eig(anal_eigvals, -1.0/(h*h), 2.0/(h*h), n);

    // finding eigenvalues with armadillo
    eig_sym(comp_eigvals, A);

    anal_eigvals = sort(anal_eigvals);
    comp_eigvals = sort(comp_eigvals);

    // testing results
    for (int i = 0; i < n; i++) {
      REQUIRE(comp_eigvals[i] == Approx(anal_eigvals[i]));
    }
}
