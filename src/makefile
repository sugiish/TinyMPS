CXX = g++
CXXFLAGS = -O2 -Wall


.PHONY: clean

mps: main.o particles.o
	g++ -Wall main.o particles.o

main.o: main.cpp
	g++ -Wall -c main.cpp

particles.o: particles.cpp
	g++ -Wall -c particles.cpp

particles.o: particles.h

clean:
	rm -f *.o *.out