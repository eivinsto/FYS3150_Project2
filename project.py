from subprocess import call, run
import os

wd = os.getcwd() + "/src"
print(wd)


def main():
    """
    Main function for running program.
    """
    build_cpp()
    test_cpp()


def build_cpp():
    call(["make", "all"], cwd=wd)


def test_cpp():
    call(["make", "test"], cwd=wd)
    run("./test_main.exe", cwd=wd)


def clean():
    call(["make", "clean"], cwd=wd)


if __name__ == '__main__':
    main()
