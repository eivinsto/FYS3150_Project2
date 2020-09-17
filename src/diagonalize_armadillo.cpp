#include <armadillo>
#include <iostream>
#include "functions.h"

using namespace std;
using namespace arma;

vec diagonalize_arma(mat A) {
  vec eigval = eig_sym(A);
  return eigval;
}
