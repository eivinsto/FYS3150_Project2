from subprocess import run
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

wd = os.getcwd() + "/src"


def build_cpp():
    run(["make", "all"], cwd=wd)


def test_cpp():
    run(["make", "test"], cwd=wd)
    run("./test_main.exe", cwd=wd)


def clean():
    run(["make", "cleandat"], cwd=wd)


choose_run = input("Choose test / toeplitz / single / double: ")
N = int(input("Size of matrix N = "))

if choose_run == "test":
    test_cpp()

if choose_run == "toeplitz":
    rho_min = float(input("rho_min = "))
    rho_max = float(input("rho_max = "))

    build_cpp()
    run(
        ["./main.exe", f"{N}", choose_run, f"{rho_min:f}", f"{rho_max:f}"],
        cwd=wd
    )

    A = np.genfromtxt(wd+f"/comp_eigvals_{N}.dat", skip_header=2)
    R = np.genfromtxt(wd+f"/comp_eigvecs_{N}.dat", skip_header=2)
    comp_eigvals = A.diagonal()
    comp_inx = comp_eigvals.argsort(kind="stable")

    anal_eigvals = np.genfromtxt(wd+f"/anal_eigvals_{N}.dat", skip_header=2)
    anal_eigvecs = np.genfromtxt(wd+f"/anal_eigvecs_{N}.dat", skip_header=2)
    clean()

    eigval = comp_eigvals[comp_inx[0]]
    rho = np.linspace(rho_min, rho_max, N-1)
    eigvec = R[:, comp_inx[0]]

    mpl.rcParams.update({"text.usetex": True})  # using latex.
    comp_lab = f"Jacobi result for {N}x{N} matrix."
    an_lab = r"Analytic result $\lambda_{1} = $ " + f"{anal_eigvals[0]}"

    plt.figure()
    plt.plot(rho, eigvec, label=comp_lab+r" $\lambda_{1} = $ " + f"{eigval}")
    plt.plot(rho, anal_eigvecs[:, 0], label=an_lab)
    plt.title(f"Eigenvector of smalest eigenvalue of {N}x{N} Toeplitz-matrix.")
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"Eigenvector $u_{1}$ as function of $\rho$.")
    plt.legend()
    plt.grid()
    plt.show()

if choose_run == "single":
    N = int(input("Size of matrix N = "))
    rho_max = float(input("rho_max = "))

    build_cpp()
    run(
        ["./main.exe", f"{N}", choose_run, f"{rho_max:f}"],
        cwd=wd
    )

    comp_eigvals = np.genfromtxt(wd+f"/quantum_eigvals_{N}.dat", skip_header=2)
    R = np.genfromtxt(wd+f"/quantum_eigvecs_{N}.dat", skip_header=2)
    clean()

    comp_inx = comp_eigvals.argsort(kind="stable")
    print(comp_eigvals[comp_inx[:4]])
