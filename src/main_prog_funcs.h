// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H

#include "sph.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace po = boost::program_options;

SPH initialise(int &n, int &T, double &h, double &dt);
void init_output_files(std::ofstream &vOut, std::ofstream &vOut2);
void time_integration(SPH &sph, int n, int T, double h, double dt,
                      std::ofstream &vOut, std::ofstream &vOut2);

#endif