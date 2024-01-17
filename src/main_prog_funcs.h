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

SPH initialise(int &n, int &T, double &h, double &dt, double &gas_constant,
               double &density_resting, double &viscosity,
               double &acceleration_gravity, double &coeff_restitution,
               double &left_wall, double &right_wall, double &bottom_wall,
               double &top_wall, int &frequency);
std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string folderPath);
void time_integration(SPH &sph, int n, int T, double h, double dt,
                      int frequency, std::ofstream &vOut, std::ofstream &vOut2);
void createDirectory(std::string folderPath);
void storeToFile(SPH &sph, int nb_particles, std::string type,
                 std::ofstream &targetFile, double dt = 0.0,
                 int currentIteration = 0);
#endif