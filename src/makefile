CXX = g++
CXXFLAGS = -O2 -Wall


.PHONY: clean

mps: main.o particles.o timer.o
	g++ -Wall main.o particles.o timer.o

main.o: main.cpp
	g++ -Wall -c main.cpp

particles.o: particles.cpp
	g++ -Wall -c particles.cpp

particles.o: particles.h

timer.o: timer.cpp
	g++ -Wall -c timer.cpp

grid.o: grid.cpp
	g++ -Wall -c grid.cpp

clean:
	rm -f *.o *.out