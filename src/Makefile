CXXFLAGS =-std=c++17

LDLIBS = -lboost_program_options -Wall -O3

default: SPH-SOLVER

SPH-main.o: SPH-main.cpp class.h initial_conditions.h
	$(CXX) $(CXXFLAGS) -c SPH-main.cpp -o SPH-main.o

class.o: class.cpp
	$(CXX) $(CXXFLAGS) -c class.cpp -o class.o	

initial_conditions.o: initial_conditions.cpp
	$(CXX) $(CXXFLAGS) -c initial_conditions.cpp -o initial_conditions.o	

SPH-SOLVER: SPH-main.o class.o initial_conditions.o	
	$(CXX) $(CXXFLAGS) -o SPH-SOLVER *.o $(LDLIBS)

.PHONY: clean

clean:
	rm SPH-SOLVER Energy-File.txt Positions-x-y.txt class.o SPH-main.o