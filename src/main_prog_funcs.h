// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H

#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

#include "fluid.h"
#include "sph_solver.h"

namespace po = boost::program_options;

fluid initialise(sph_solver &solver);
std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string folderPath);

void createDirectory(std::string folderPath);
void storeToFile(fluid &fluid, std::string type,
                 std::ofstream &targetFile, double dt = 0.0,
                 int currentIteration = 0);
#endif