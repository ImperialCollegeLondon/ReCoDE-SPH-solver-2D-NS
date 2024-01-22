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

struct SimulationParameters {
  // Declare the parameters of the problem
  int total_iter;  // total number of iterations required for the time
                   // integration
  double dt;       // time step

  // Constants
  double h;
  double gas_constant;
  double density_resting;
  double viscosity;
  double acceleration_gravity;
  double coeff_restitution;

  // Domain boundaries
  double left_wall;
  double right_wall;
  double bottom_wall;
  double top_wall;

  // Number of particles
  int nb_particles;

  // Output frequency
  int frequency;
};

SPH initialise(SimulationParameters &parameters);
std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string folderPath);
void time_integration(SPH &sph, SimulationParameters simulationParameters,
                      std::ofstream &vOut, std::ofstream &vOut2);
void createDirectory(std::string folderPath);
void storeToFile(SPH &sph, int nb_particles, std::string type,
                 std::ofstream &targetFile, double dt = 0.0,
                 int currentIteration = 0);
#endif