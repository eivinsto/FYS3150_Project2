CC := g++
FLAGS := -Wall -O2 -march=native
LIBS := -larmadillo
OFILES := main.o

.PHONY: all clean run

all: main.exe

main.exe: $(OFILES)
	$(CC) $(FLAGS) $(OFILES) -o main.exe $(LIBS)

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp -c

clean:
	rm -f *.o
	rm -f *.exe
	rm -f *.bin

cleanbin:
	rm -f *.bin

run:
	./main.exe

print-% : ; @echo $* = $($*)
