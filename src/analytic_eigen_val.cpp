#include <armadillo>
#include <cmath>
#include "functions.h"

using namespace arma;
using namespace std;

vec anal_eig(double d, double a, int N) {
  vec eigval = zeros<vec>(N);
  double angval1 = 2*M_PI/N;

  for (int i = 0; i < N; i++) {
    double angval2 = (i + 1)*angval1;
    eigval[i] = 2*a*cos(angval2);
  }

  return eigval;
}
