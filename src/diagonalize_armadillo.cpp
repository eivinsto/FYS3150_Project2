#include <armadillo>
#include <iostream>
#include "functions.h"

using namespace std;
using namespace arma;

void diagonalize_arma(vec& eigval, mat A) {
  eig_sym(eigval, A);
}
