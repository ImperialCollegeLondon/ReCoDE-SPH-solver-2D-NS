// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H

#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>

#include "sph.h"

namespace po = boost::program_options;

SPH initialise(int &n, int &T, double &h, double &dt);
std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string folderPath);
void time_integration(SPH &sph, int n, int T, double h, double dt, int freq,
                      std::ofstream &vOut, std::ofstream &vOut2);
void createDirectory(std::string folderPath);
void storeToFile(SPH &sph, int nb_particles, std::string type,
                 std::ofstream &targetFile, double dt = 0.0,
                 int currentIteration = 0);
#endif