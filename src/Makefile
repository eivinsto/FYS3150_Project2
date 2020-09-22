CC := g++
FLAGS := -Wall -O2 -march=native
LIBS := -larmadillo
OFILES := main.o create_tridiagonal_matrix.o diagonalize_armadillo.o jacobi_solver.o analytic_eigen_val.o

.PHONY: all clean run

all: main.exe

main.exe: $(OFILES)
	$(CC) $(FLAGS) $(OFILES) -o main.exe $(LIBS)

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp -c

create_tridiagonal_matrix.o: create_tridiagonal_matrix.cpp
	$(CC) $(FLAGS) create_tridiagonal_matrix.cpp -c

diagonalize_armadillo.o: diagonalize_armadillo.cpp
	$(CC) $(FLAGS) diagonalize_armadillo.cpp -c

jacobi_solver.o: jacobi_solver.cpp
	$(CC) $(FLAGS) jacobi_solver.cpp -c

analytic_eigen_val.o: analytic_eigen_val.cpp
	$(CC) $(FLAGS) analytic_eigen_val.cpp -c

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.bin

cleanbin:
	rm -f *.bin

run:
	./main.exe

print-% : ; @echo $* = $($*)
