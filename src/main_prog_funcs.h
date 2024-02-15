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
  double total_time;  // total time of the simulation
  int total_iter;     // total number of iterations required for the time
                      // integration
  double dt;          // time step

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

// Fixed nb_particles ic cases
std::map<std::string, int> initConditionToParticlesMap = {
    {"ic-one-particle", 1},
    {"ic-two-particles", 2},
    {"ic-three-particles", 3},
    {"ic-four-particles", 4}};

SPH initialise(SimulationParameters &parameters);
std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string folderPath);
void time_integration(SPH &sph,
                      const SimulationParameters &simulationParameters,
                      std::ofstream &vOut, std::ofstream &vOut2);
void retrieveInputsFromFile(std::string fileName, std::string icCase,
                            po::options_description desc,
                            po::variables_map &vm);
void handleInputErrors(SimulationParameters input);
SPH setInitialConditions(std::string ic_case,
                         SimulationParameters &simulationParameters,
                         po::variables_map ic_vm);
void createDirectory(std::string folderPath);
void storeToFile(SPH &sph, int nb_particles, std::string type,
                 std::ofstream &targetFile, double dt = 0.0,
                 int currentIteration = 0);
#endif