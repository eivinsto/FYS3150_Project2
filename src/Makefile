CC := g++
FLAGS := -Wall -O2 -march=native
LIBS := -larmadillo
OFILES := jacobi_solver.o functions.o

.PHONY: all clean cleanbin test

all: main.exe

main.exe: main.o $(OFILES)
	$(CC) $(FLAGS) main.o $(OFILES) -o main.exe $(LIBS)


main.o: main.cpp
	$(CC) $(FLAGS) main.cpp -c

jacobi_solver.o: jacobi_solver.cpp
	$(CC) $(FLAGS) jacobi_solver.cpp -c

functions.o: functions.cpp
	$(CC) $(FLAGS) functions.cpp -c

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.bin

cleanbin:
	rm -f *.bin

test: test_main.exe

test_main.exe: test_main.o $(OFILES) test_functions.o
	$(CC) $(FLAGS) test_main.o $(OFILES) test_functions.o -o test_main.exe $(LIBS)

test_main.o: test_main.cpp
	$(CC) $(FLAGS) test_main.cpp -c

test_functions.o: test_functions.cpp
	$(CC) $(FLAGS) test_functions.cpp -c

print-% : ; @echo $* = $($*)
