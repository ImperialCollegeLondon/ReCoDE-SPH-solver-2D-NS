// Definition of the functions called by the main program
#ifndef MAIN_H
#define MAIN_H

#include "sph_2d.h"
#include <boost/program_options.hpp>
#include <iostream>

namespace po = boost::program_options;

sph_2d initialise(int &n, int &T, double &h, double &dt);
void init_output_files(std::ofstream &vOut, std::ofstream &vOut2);

#endif