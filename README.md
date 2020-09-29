# FYS3150_Project2
All code for this report was written in python 3.8, and C++.
To generate the data used in the report, run the python-script "project.py" in the root directory of this repo.

## Usage:
When executed, this script will ask you which part of the report to run:
 * te / test     - Run unit-tests for the jacobi_solver class.
 * b / benchmark - Runs benchmark of the jacobi_solver and Armadillo eigen_sym solver for comparison.
 * to / toeplitz - Solves the buckling beam problem, and plots the analytic and numerical eigenvector for the smallest eigenvalue.
 * s / single    - Solves 4 lowest eigenvalues of single-electron atom, and compares them with analytic values.
 * d / double    - Finds ground-state energies of a two-electron atom, and compares them with analytic approximate values.

Example run:
```console
$ python project.py
To select calculation, write whats in brackets, or whole word from list of options below:
[te]st / [b]enchmark / [to]eplitz / [s]ingle / [d]ouble
Choose run: d
Size of matrix N = 400
rho_max = 115.625
make: Nothing to be done for 'all'.
rm -f *.dat
Eigenvalues of two-electron atom. N = 400, rho_max = 115.625
Omega_r:    Analytic:    Numerical:    Relative error:
    0.01            0         0.106          6.966e-03
    0.50            2         2.189          6.434e-02
    1.00            4         3.898          7.623e-02
    5.00           14        14.118          4.817e-03
```

When prompted to enter value for N, you give the N used in the report
methodology. Hence, N = 100 will create matrices of size (N-1)x(N-1) -> 100x100.
For the buckling beam problem, the plotted eigenvectors will have length N+1, with the first and last element fixed to zero.
