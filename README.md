# FYS3150_Project2
All code for this report was written in python 3.8, and C++.
To generate the data used in the report, run the python-script "project.py" in the root directory of this repo.

## Usage:
When executed, this script will ask you which part of the report to run:
 * te / test     - Run unit-tests for the jacobi\_solver class.
 * b / benchmark - Runs benchmark of the jacobi\_solver and Armadillo eigen\_sym solver for comparison.
 * to / toeplitz - Solves the buckling beam problem, and plots the analytic and numerical eigenvector for the smallest eigenvalue.
 * s / single    - Solves 4 lowest eigenvalues of single-electron atom, and compares them with analytic values.
 * d / double    - Finds ground-state energies of a two-electron atom, and compares them with analytic approximate values.

Example run:
```console
$ python project.py
To select calculation, write whats in brackets, or whole word from list of options below:
[te]st / [b]enchmark / [to]eplitz / [s]ingle / [d]ouble
Choose run: s
Size of matrix N = 400
rho_max = 25
make: Nothing to be done for 'all'.
rm -f *.dat
Eigenvalues of single-electron atom. N = 400, rho_max = 25.0
Analytic:    Numerical:    Relative error:
       3        2.999          4.071e-04
       7        6.994          8.727e-04
      11       10.985          1.356e-03
      15       14.972          1.843e-03
```

When prompted to enter value for N, you give the N used in the report
methodology. Hence, N = 100 will create matrices of size (N-1)x(N-1) = 99x99.
For the buckling beam problem, the plotted eigenvectors will have length N+1, with the first and last element fixed to zero.
