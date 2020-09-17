#include <armadillo>
#include <iostream>
#include "functions.h"

using namespace std;
using namespace arma;

cx_vec diagonalize_arma(mat A) {
  cx_vec eigval;
  eig_gen(eigval, A, "balance");

  return eigval;
}
