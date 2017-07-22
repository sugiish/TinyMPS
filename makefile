CXX = g++
CXXFLAGS = -O2 -Wall

COMPILE.C = $(COMPILE.cc)
COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c

OUTPUT_OPTION = -o $@
SOURCE_DIR = src
OBJECT_DIR = obj

VPATH = src include
CPPFLAGS = -I include

mps.out: main.o particles.o timer.o
	g++ -Wall $(CPPFLAGS) $^ -o $@

main.o: main.cpp
	g++ -Wall $(CPPFLAGS) -c $<

particles.o: particles.cpp
	g++ -Wall $(CPPFLAGS) -c $<

particles.o: particles.h

timer.o: timer.cpp
	g++ -Wall $(CPPFLAGS) -c $<

grid.o: grid.cpp
	g++ -Wall $(CPPFLAGS) -c $<

.PHONY: clean

clean:
	rm -f *.o *.out

#%.o: %.c
#	$(COMPILE.CC) $(OUTPUT_OPTION) $<