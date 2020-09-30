"""
Python script to generate the data used in the report.

When executed, this script will ask you which part of the report to run:
    te / test     - Run unit-tests for the jacobi_solver class.
    b / benchmark - Runs benchmark of the jacobi_solver and Armadillo
                    eigen_sym solver for comparison.
    to / toeplitz - Solves the buckling beam problem, and plots the
                    analytic and numerical eigenvector for the smallest
                    eigenvalue.
    s / single    - Solves 4 lowest eigenvalues of single-electron quantum dot,
                    and compares them with analytic values.
    d / double    - Finds ground-state energies of a two-electron quantum dot, and
                    compares them with analytic approximate values.

When prompted to enter value for N, you give the N used in the report
methodology. Hence, N = 100 will create matrices of size (N-1)x(N-1)
as outlined in the report.
"""
from subprocess import run
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import shutil

# retriveing working directory:
pwd = os.getcwd()
wd = pwd + "/src"


def build_cpp():
    """Function for building c++ program."""
    run(["make", "all"], cwd=wd)


def test_cpp():
    """Function for running unit-tests."""
    run(["make", "test"], cwd=wd)
    run("./test_main.exe", cwd=wd)


def clean():
    """Function for cleaning datafiles in src directory."""
    run(["make", "cleandat"], cwd=wd)


print("To select calculation, write whats in brackets", end="")
print(", or whole word from list of options below:")
print("[te]st / [b]enchmark / [to]eplitz / [s]ingle / [d]ouble")
choose_run = input("Choose run: ")  # asking user for what to do.


if choose_run == "test" or choose_run == "te":
    """Performing unit-tests."""
    test_cpp()


if choose_run == "toeplitz" or choose_run == "to":
    """Solving buckling beam problem."""
    N = int(input("Size of matrix N = "))  # retriveing N.
    rho_min = float(input("rho_min = "))  # retriveing max rho value.
    rho_max = float(input("rho_max = "))  # retriveing min rho value.

    # building and running c++ program:
    build_cpp()
    run(
        ["./main.exe", f"{N}", choose_run, f"{rho_min:f}", f"{rho_max:f}"],
        cwd=wd
    )

    # retriveing numerical solution from files:
    A = np.genfromtxt(wd+f"/comp_eigvals_{N}.dat", skip_header=2)
    R = np.genfromtxt(wd+f"/comp_eigvecs_{N}.dat", skip_header=2)
    comp_eigvals = A.diagonal()  # getting eigenvalues from diagonal matrix.
    comp_inx = comp_eigvals.argsort(kind="stable")  # sorting eigenvalues.

    # retriveing analytic values from files:
    anal_eigvals = np.genfromtxt(wd+f"/anal_eigvals_{N}.dat", skip_header=2)
    anal_eigvecs = np.genfromtxt(wd+f"/anal_eigvecs_{N}.dat", skip_header=2)

    # copying data to data directory:
    shutil.move(wd+f"/anal_eigvals_{N}.dat",
                pwd + "/data/toeplitz_anal_eigvals.dat")
    shutil.move(wd+f"/anal_eigvecs_{N}.dat",
                pwd + "/data/toeplitz_anal_eigvecs.dat")
    shutil.move(wd+f"/comp_eigvals_{N}.dat",
                pwd + "/data/toeplitz_comp_eigvals.dat")
    shutil.move(wd+f"/comp_eigvecs_{N}.dat",
                pwd + "/data/toeplitz_comp_eigvecs.dat")
    clean()  # cleaning data files from src directory.

    eigval = comp_eigvals[comp_inx[0]]  # retriveing smallest eigenvalue.
    rho = np.linspace(rho_min, rho_max, N+1)  # creating rho values.
    eigvec = np.zeros(N+1)  # setting boundary conditions.
    eigvec[1:-1] = R[:, comp_inx[0]]  # retriveing solution.

    anal_eigvec = np.zeros(N+1)  # setting boundary conditions.
    anal_eigvec[1:-1] = anal_eigvecs[:, 0]  # retriveing analytic solution.

    # calculating the norms of the eigenvectors:
    anal_norm = np.linalg.norm(anal_eigvec)
    comp_norm = np.linalg.norm(eigvec)

    mpl.rcParams.update({"text.usetex": True})  # using latex.
    # creating labels for plot:
    comp_lab = r"Jacobi result $\lambda_{1} = $ " + f"{eigval:.3f}"
    an_lab = r"Analytic result $\lambda_{1} = $ " + f"{anal_eigvals[0]:.3f}"

    # plotting normalized numerical and analytic eigenvectors:
    plt.figure()
    plt.plot(rho, eigvec/comp_norm, label=comp_lab)
    plt.plot(rho, anal_eigvec/anal_norm, label=an_lab)
    plt.title(f"Eigenvector of smallest eigenvalue of Toeplitz-matrix. {N = }")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"Eigenvector $\vec{u}_{1}$ as function of $\rho$.")
    plt.legend()
    plt.grid()
    plt.show()


if choose_run == "single" or choose_run == "s":
    """Finding 4 lowest energies of one-electron system."""
    anal_eigvals = np.array([3, 7, 11, 15])  # analytic solution.

    N = int(input("Size of matrix N = "))  # retriveing N.
    rho_max = float(input("rho_max = "))  # retriveing rho_max.

    # building and running c++ program:
    build_cpp()
    run(
        ["./main.exe", f"{N}", choose_run, f"{rho_max:f}"],
        cwd=wd
    )

    # retriveing numerical results from file:
    comp_eigvals = np.genfromtxt(wd+f"/quantum_eigvals_{N}.dat", skip_header=2)
    clean()  # cleaning files in src directory.

    #
    comp_eigvals.sort(kind="stable")  # sorting eigenvalues.

    # calculating relative error:
    err = np.abs(anal_eigvals-comp_eigvals[:4])/anal_eigvals

    # writing results to file in data directory:
    with open(pwd + "/data/single_electron_data.dat", "w") as output:
        header1 = f"Eigenvalues of single-electron quantum dot. {N = }, {rho_max = }"
        header2 = "Analytic:    Numerical:    Relative error:"
        print(header1)
        output.write(header1 + "\n")

        print(header2)
        output.write(header2 + "\n")

        for i in range(len(err)):
            analstr = f"{anal_eigvals[i]:8.0f}"
            numstr = f"{comp_eigvals[i]:9.3f}"
            errstr = f"{err[i]:15.3e}"
            line = analstr + "    " + numstr + "    " + errstr
            print(line)
            output.write(line + "\n")


if choose_run == "double" or choose_run == "d":
    """
    Finds ground-state energies of two-electron quantum dot with
    different values of omega_r.
    """
    N = int(input("Size of matrix N = "))  # retriveing N.
    rho_max = float(input("rho_max = "))  # retriveing rho_max.

    omega_r = np.array([0.01, 0.5, 1, 5])  # omega_r values.
    N_omega = len(omega_r)  # number of omega_r values.

    # empty arrays for numerical and analytic ground-state energies:
    num_eigvals = np.empty((N-1, N_omega))
    anal_eigvals = np.empty((N-1, N_omega))
    err = np.empty(N_omega)  # empty array for relative errors.

    build_cpp()  # building c++ program.

    for i in range(N_omega):
        # running c++ program for each omega_r value:
        run(
            [
                "./main.exe",
                f"{N}",
                choose_run,
                f"{rho_max:f}",
                f"{omega_r[i]:f}"
            ],
            cwd=wd
        )

        # retriveing numerical and analytic energies for each omega_r:
        anal_eigvals[:, i] = np.genfromtxt(
            wd + f"/double_qdot_anal_eigvals_{N}.dat", skip_header=2
        )

        num_eigvals[:, i] = np.genfromtxt(
            wd + f"/double_qdot_num_eigvals_{N}.dat", skip_header=2
        )

        num_eigvals[:, i] = np.sort(num_eigvals[:, i], kind="stable")

        # calculating relative error for each omega_r:
        err[i] = np.abs(
            (anal_eigvals[0, i]-num_eigvals[0, i]) /
            anal_eigvals[0, i]
        )

    clean()  # cleaning data files from src directory.

    # writing results to file in data directory:
    with open(pwd + "/data/double_electron_data.dat", "w") as output:
        header1 = f"Eigenvalues of two-electron quantum dot. {N = }, {rho_max = }"
        header2 = "Omega_r:    Analytic:    Numerical:    Relative error:"
        print(header1)
        output.write(header1 + "\n")

        print(header2)
        output.write(header2 + "\n")

        for i in range(len(omega_r)):
            omegastr = f"{omega_r[i]:8.2f}"
            analstr = f"{anal_eigvals[0, i]:9.3f}"
            numstr = f"{num_eigvals[0, i]:10.3f}"
            errstr = f"{err[i]:15.3e}"

            line = (
                omegastr + "    " +
                analstr + "    " +
                numstr + "    " +
                errstr
            )

            print(line)
            output.write(line + "\n")


if choose_run == "benchmark" or choose_run == "b":
    """
    Performing benchmark of jacobi method,
    with Armadillo eig_sym for comparison.
    """
    N = 501

    # builds and runs c++ program:
    build_cpp()
    print()
    run(["./main.exe", f"{N}", choose_run], cwd=wd)
    print()

    # retriveing benchmark results:
    benchmark_times = np.genfromtxt(wd + "/benchmark_times.dat", skip_header=2)
    clean()  # clean data from src directory.

    # calculating mean and standard deviation for Armadillo:
    arma_toeplitz_mean = np.mean(benchmark_times[:, 0])
    arma_toeplitz_std = np.std(benchmark_times[:, 0])
    arma_qdot_mean = np.mean(benchmark_times[:, 2])
    arma_qdot_std = np.std(benchmark_times[:, 2])

    # calculating mean and standard deviation for jacobi_solver:
    jacobi_toeplitz_mean = np.mean(benchmark_times[:, 1])
    jacobi_toeplitz_std = np.std(benchmark_times[:, 1])
    jacobi_qdot_mean = np.mean(benchmark_times[:, 3])
    jacobi_qdot_std = np.std(benchmark_times[:, 3])

    # printing results:
    header1 = f"Time spent solving {N-1}x{N-1} Toeplitz matrix."
    header2 = f"{N = }, rho_max = 1, rho_min = 0"
    jacobistr = f"Jacobi solver: {jacobi_toeplitz_mean:.4f} s \u00B1 {jacobi_toeplitz_std:.4f} s"
    armastr = f"Armadillo solver: {arma_toeplitz_mean:.4f} s \u00B1 {arma_toeplitz_std:.4f} s"
    print(header1)
    print(header2)
    print(jacobistr)
    print(armastr)

    # writing results to file in data directory:
    with open(pwd + "/data/benchmark_toeplitz.dat", "w") as output:
        output.write(header1 + "\n")
        output.write(header2 + "\n")
        output.write(jacobistr + "\n")
        output.write(armastr + "\n")

    # printing results:
    header1 = f"Time spent solving {N-1}x{N-1} two-electron matrix."
    header2 = f"{N = }, rho_max = 144.53125, omega_r = 5.0"
    jacobistr = f"Jacobi solver: {jacobi_qdot_mean:.4f} s \u00B1 {jacobi_qdot_std:.4f} s"
    armastr = f"Armadillo solver: {arma_qdot_mean:.4f} s \u00B1 {arma_qdot_std:.4f} s"
    print(header1)
    print(header2)
    print(jacobistr)
    print(armastr)

    # writing results to file in data directory:
    with open(pwd + "/data/benchmark_double.dat", "w") as output:
        output.write(header1 + "\n")
        output.write(header2 + "\n")
        output.write(jacobistr + "\n")
        output.write(armastr + "\n")
