// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H
#include <memory>

#include "fluid.h"
#include "sph_solver.h"

void initialise(std::unique_ptr<Fluid> &fluidPtr, SphSolver &solver);
std::tuple<std::ofstream, std::ofstream, std::ofstream> initOutputFiles(
    std::string folderPath);

void createDirectory(std::string folderPath);

void storeToFile(Fluid &fluid, std::string type, std::ofstream &targetFile,
                 double dt = 0.0, int currentIteration = 0);
#endif