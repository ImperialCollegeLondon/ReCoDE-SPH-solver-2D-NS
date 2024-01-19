// This is a main program to run the SPH simulation
#include <boost/program_options.hpp>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <tuple>

#include "initial_conditions.h"
#include "main_prog_funcs.h"
#include "sph.h"

// Start of the main programme
int main(int argc, char *argv[]) {
  // Declare the parameters of the problem
  int nb_particles;  // number of particles
  int total_iter;    // total number of iterations required for the time
                     // integration
  double h;
  double dt;
  std::string OUTPUT_FOLDER = "../output";
  // Read input files, initialise the sph class and the parameters of the
  // problem
  SPH sph = initialise(nb_particles, total_iter, h, dt);

  std ::cout << "Initialisation finished -- OK"
             << "\n";

  // Initialise output files
  auto [initialPositions, finalPositionsFile, energiesFile] =
      init_output_files(OUTPUT_FOLDER);

  std ::cout << "Output files created -- OK"
             << "\n";

  // Store particles' positions before integration has started
  storeToFile(sph, nb_particles, "position", initialPositions);
  initialPositions.close();

  int frequency = 3;  //! Delete this when frequency is read from input variable
  // Time integration loop
  time_integration(sph, nb_particles, total_iter, h, dt, frequency,
                   finalPositionsFile, energiesFile);

  std ::cout << "SPH-SOLVER executed successfully -- OK"
             << "\n";

  return 0;
}

SPH initialise(int &nb_particles, int &total_iter, double &h, double &dt) {
  // Process to obtain the directions provided by the user
  po::options_description desc("Allowed options");
  desc.add_options()("init_condition", po::value<std::string>(),
                     "take an initial condition")("T", po::value<double>(),
                                                  "take integration time")(
      "dt", po::value<double>(), "take time-step")("h", po::value<double>(),
                                                   "take radius of influence");

  po::variables_map vm;
  std::ifstream inputFile;
  inputFile.open("../input/case.txt");

  if (inputFile.is_open()) {
    po::store(po::parse_config_file(inputFile, desc), vm);
    inputFile.close();
  } else {
    std::cerr << "Error opening file: case.txt" << std::endl;
  }

  po::notify(vm);

  int n1 = 17;   // required for ic-block-drop
  int n2 = 25;   // required for ic-block-drop
  int n3 = 100;  // required for ic-dam-break and ic-droplet

  double total_time = vm["T"].as<double>();  // Total integration time
  dt = vm["dt"].as<double>();                // Time step dt
  h = vm["h"].as<double>();                  // Radius of influence

  total_iter =
      ceil(total_time / dt);  // Transform time in seconds to iterations

  std::map<std::string, int> initConditionToParticlesMap = {
      {"ic-one-particle", 1},      {"ic-two-particles", 2},
      {"ic-three-particles", 3},   {"ic-four-particles", 4},
      {"ic-dam-break", n3},        {"ic-block-drop", n1 * n2},
      {"ic-droplet", dropletn(n3)}};

  nb_particles =
      initConditionToParticlesMap[vm["init_condition"].as<std::string>()];

  // Define the solver object (called sph)
  // In its definition, the number of particles is required
  SPH sph(nb_particles);

  /**After the number of particles is introduced inside the class and
   * therefore the appropriate matrices are initialized, the particles
   * are ordered in the correct positions
   **/

  // Create map to associate function names with function pointers
  std::map<std::string, std::function<void(int, SPH &)>> functionMap = {
      {"ic-one-particle", ic_one_particle},
      {"ic-two-particles", ic_two_particles},
      {"ic-three-particles", ic_three_particles},
      {"ic-four-particles", ic_four_particles},
      {"ic-dam-break", ic_dam_break},
      {"ic-droplet", ic_droplet}};

  // Get the function pointer from the map
  auto initFunc = functionMap.find(vm["init_condition"].as<std::string>());

  if (initFunc != functionMap.end()) {
    int n_particles = nb_particles;

    // The ic-droplet case requires a different n argument.
    if (vm["init_condition"].as<std::string>() == "ic-droplet") {
      n_particles = n3;
    }

    // Retrieves and runs the provided function object
    initFunc->second(n_particles, sph);

  } else if (vm["init_condition"].as<std::string>() == "ic-block-drop") {
    /**The ic-block-drop case is not in the map because it has two
     * additional parameters, so it requires a different case.
     **/
    ic_block_drop(nb_particles, n1, n2, sph);
  } else {
    std::cerr << "Error: Function not found!" << std::endl;
  }

  sph.set_timestep(dt);
  sph.set_rad_infl(h);

  // Calculate the mass of the particles
  sph.calc_mass();

  return sph;
}

void createDirectory(std::string folderPath) {
  // Check if the target folder already exists
  if (!std::filesystem::exists(folderPath)) {
    // Create target folder
    std::filesystem::create_directories(folderPath);
  }
}

std::tuple<std::ofstream, std::ofstream, std::ofstream> init_output_files(
    std::string outputFolder) {
  // Create the output folder if it doesn't exist
  createDirectory(outputFolder);

  // Declare and initialise the output files
  std::ofstream initialPositions(outputFolder + "/initial-positions.csv",
                                 std::ios::out | std::ios::trunc);
  std::ofstream finalPositions(outputFolder + "/final-positions.csv",
                               std::ios::out | std::ios::trunc);
  std::ofstream energies(outputFolder + "/energies.csv",
                         std::ios::out | std::ios::trunc);

  initialPositions << std::fixed << std::setprecision(5);
  initialPositions << "Position_X,Position_Y"
                   << "\n";

  finalPositions << std::fixed << std::setprecision(5);
  finalPositions << "Position_X,Position_Y"
                 << "\n";

  energies << std::fixed << std::setprecision(5);
  energies << "t,Ek,Ep,Etotal"
           << "\n";

  return std::make_tuple(std::move(initialPositions), std::move(finalPositions),
                         std::move(energies));
}

void storeToFile(const SPH &sph, int nb_particles, const std::string &type,
                 std::ofstream &targetFile, double dt, int currentIteration) {
  if (type == "energy") {
    // Write energies on the Energy-File
    targetFile << currentIteration * dt << "," << sph.return_kinetic_energy()
               << "," << sph.return_potential_energy() << ","
               << sph.return_potential_energy() + sph.return_kinetic_energy()
               << "\n";
  } else if (type == "position") {
    for (int k = 0; k < nb_particles; k++) {
      targetFile << sph.return_position_x(k) << "," << sph.return_position_y(k)
                 << "\n";
    }
  }
}

void time_integration(SPH &sph, int nb_particles, int total_iter, double h,
                      double dt, int frequency,
                      std::ofstream &finalPositionsFile,
                      std::ofstream &energiesFile) {
  std ::cout << "Time integration started -- OK"
             << "\n";

  for (int t = 0; t < total_iter; t++) {
    // In each iteration the distances between the particles are recalculated,
    // as well as their densities
    sph.calc_particle_distance();
    sph.calc_density();
    sph.calc_pressure();
    sph.particle_iterations();

    if (t % frequency == 0) {
      storeToFile(sph, nb_particles, "energy", energiesFile, dt, t);
    }
  }
  // Store particles' positions after integration is completed
  storeToFile(sph, nb_particles, "position", finalPositionsFile);

  std ::cout << "Time integration finished -- OK"
             << "\n";
}
