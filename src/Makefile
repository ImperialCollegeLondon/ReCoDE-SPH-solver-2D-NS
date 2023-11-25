CXX=mpicxx
CXXFLAGS =-std=c++11

LDLIBS = -lboost_program_options -Wall -O3

default: SPH-SOLVER

SPH-main.o: SPH-main.cpp class.h
	$(CXX) $(CXXFLAGS) -c SPH-main.cpp -o SPH-main.o

class.o: class.cpp
	$(CXX) $(CXXFLAGS) -c class.cpp -o class.o	

SPH-SOLVER: SPH-main.o class.o
	$(CXX) $(CXXFLAGS) -o SPH-SOLVER *.o $(LDLIBS)

.PHONY: clean

clean:
	rm SPH-SOLVER Energy-File.txt Positions-x-y.txt class.o SPH-main.o
