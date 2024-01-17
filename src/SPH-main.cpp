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
  int total_iter;  // total number of iterations required for the time
                   // integration
  double dt;

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

  int frequency;

  std::string OUTPUT_FOLDER = "../output";
  // Read input files, initialise the sph class and the parameters of the
  // problem
  SPH sph =
      initialise(nb_particles, total_iter, h, dt, gas_constant, density_resting,
                 viscosity, acceleration_gravity, coeff_restitution, left_wall,
                 right_wall, bottom_wall, top_wall, frequency);

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

  // Time integration loop
  time_integration(sph, nb_particles, total_iter, h, dt, frequency,
                   finalPositionsFile, energiesFile);

  std ::cout << "SPH-SOLVER executed successfully -- OK"
             << "\n";

  return 0;
}

SPH initialise(int &nb_particles, int &total_iter, double &h, double &dt,
               double &gas_constant, double &density_resting, double &viscosity,
               double &acceleration_gravity, double &coeff_restitution,
               double &left_wall, double &right_wall, double &bottom_wall,
               double &top_wall, int &frequency) {
  // Process to obtain the inputs provided by the user
  po::options_description desc("Allowed options");
  desc.add_options()("init_condition", po::value<std::string>(),
                     "take an initial condition")("T", po::value<double>(),
                                                  "take integration time")(
      "dt", po::value<double>(), "take time-step")("h", po::value<double>(),
                                                   "take radius of influence")(
      "gas_constant", po::value<double>(), "take gas constant")(
      "density_resting", po::value<double>(), "take resting density")(
      "viscosity", po::value<double>(), "take viscosity")(
      "acceleration_gravity", po::value<double>(), "take acc due to gravity")(
      "coeff_restitution", po::value<double>(), "take coeff of restitution")(
      "left_wall", po::value<double>(), "take left wall position")(
      "right_wall", po::value<double>(), "take right wall position")(
      "bottom_wall", po::value<double>(), "take bottom wall position")(
      "top_wall", po::value<double>(), "take top wall position")(
      "length", po::value<double>(), "take length of the block")(
      "width", po::value<double>(), "take width of the block")(
      "radius", po::value<double>(), "take radius of the droplet")(
      "n", po::value<int>(), "take number of particles")(
      "center_x", po::value<double>(), "take center of the particle mass in x")(
      "center_y", po::value<double>(), "take center of the particle mass in y")(
      "init_x_1", po::value<double>(), "take x_1")(
      "init_y_1", po::value<double>(), "take y_2")(
      "init_x_2", po::value<double>(), "take x_2")(
      "init_y_2", po::value<double>(), "take y_2")(
      "init_x_3", po::value<double>(), "take x_3")(
      "init_y_3", po::value<double>(), "take y_3")(
      "init_x_4", po::value<double>(), "take x_4")(
      "init_y_4", po::value<double>(), "take y_4")(
      "output_frequency", po::value<int>(), "take frequency of output");

  // Map the inputs read from the case file to expected inputs
  po::variables_map case_vm;
  std::ifstream caseFile;
  caseFile.open("../inputs/case.txt");

  if (caseFile.is_open()) {
    po::store(po::parse_config_file(caseFile, desc), case_vm);
    caseFile.close();
  } else {
    std::cerr << "Error opening file: case.txt" << std::endl;
  }
  po::notify(case_vm);

  double total_time = case_vm["T"].as<double>();  // Total integration time
  // Error handling for the total integration time
  if (total_time <= 0) {
    std::cerr << "Error: Total integration time must be positive!" << std::endl;
    exit(1);
  }

  dt = case_vm["dt"].as<double>();  // Time step dt
  // Error handling for the time step
  if (dt <= 0 or dt > total_time) {
    std::cerr << "Error: Time step must be positive and lower than the total "
              << "integration time!" << std::endl;
    exit(1);
  }

  total_iter =
      ceil(total_time / dt);  // Transform time in seconds to iterations

  frequency = case_vm["output_frequency"].as<int>();
  // Error handling for the output frequency
  if (frequency <= 0 or frequency > total_iter) {
    std::cerr << "Error: Output frequency must be positive and lower than the "
              << "total number of iterations!" << std::endl;
    exit(1);
  }

  // Map the inputs read from the domain file to expected inputs
  po::variables_map domain_vm;
  std::ifstream domainFile;
  domainFile.open("../inputs/domain.txt");

  if (domainFile.is_open()) {
    po::store(po::parse_config_file(domainFile, desc), domain_vm);
    domainFile.close();
  } else {
    std::cerr << "Error opening file: domain.txt" << std::endl;
  }

  po::notify(domain_vm);

  left_wall = domain_vm["left_wall"].as<double>();
  right_wall = domain_vm["right_wall"].as<double>();
  bottom_wall = domain_vm["bottom_wall"].as<double>();
  top_wall = domain_vm["top_wall"].as<double>();

  // Error handling for the domain boundaries
  if (left_wall >= right_wall || bottom_wall >= top_wall) {
    std::cerr << "Error: Please adjust your domain boundaries so that "
              << "left_wall < right wall and bottom_wall < top_wall."
              << std::endl;
    exit(1);
  }

  // Map the inputs read from the initial condition file to expected inputs
  std::string ic_case = case_vm["init_condition"].as<std::string>();
  po::variables_map ic_vm;
  std::ifstream icFile;
  // Open the file of the initial condition the user has chosen
  icFile.open("../inputs/" + ic_case + ".txt");

  if (icFile.is_open()) {
    po::store(po::parse_config_file(icFile, desc), ic_vm);
    icFile.close();
  } else {
    std::cerr << "Error opening file: " << ic_case << ".txt. Make sure "
              << "that the value of the init_condition in the case.txt file is "
              << "one of the following: ic-one-particle, ic-two-particles, "
              << "ic-three-particles, ic-four-particles, ic-droplet, "
              << "ic-block-drop." << std::endl;
    exit(1);
  }
  po::notify(ic_vm);

  // Fixed nb_particles ic cases
  std::map<std::string, int> initConditionToParticlesMap = {
      {"ic-one-particle", 1},
      {"ic-two-particles", 2},
      {"ic-three-particles", 3},
      {"ic-four-particles", 4}};

  // Get the number of particles based on the ic case
  if (ic_case == "ic-droplet" || ic_case == "ic-block-drop") {
    nb_particles = ic_vm["n"].as<int>();
    // Error handling for the number of particles
    if (nb_particles <= 0) {
      std::cerr << "Error: Number of particles must be positive!" << std::endl;
      exit(1);
    }
  } else {
    nb_particles = initConditionToParticlesMap[case_vm["init_condition"]
                                                   .as<std::string>()];
  }
  /**After the number of particles is introduced inside the class and
   * therefore the appropriate matrices are initialized, the particles
   * are ordered in the correct positions
   **/
  SPH sph(nb_particles);
  if (ic_case == "ic-one-particle" || ic_case == "ic-two-particles" ||
      ic_case == "ic-three-particles" || ic_case == "ic-four-particles") {
    double *init_x = new double[nb_particles];
    double *init_y = new double[nb_particles];
    for (int i = 0; i < nb_particles; i++) {
      init_x[i] = ic_vm["init_x_" + std::to_string(i + 1)].as<double>();
      init_y[i] = ic_vm["init_y_" + std::to_string(i + 1)].as<double>();
      // Error handling for the initial positions
      if (init_x[i] < left_wall || init_x[i] > right_wall ||
          init_y[i] < bottom_wall || init_y[i] > top_wall) {
        std::cerr << "Error: Particles must be within the domain boundaries! "
                  << "Please adjust the initial position coordinates."
                  << std::endl;
        exit(1);
      }
    }
    sph = ic_basic(nb_particles, init_x, init_y);
  } else if (ic_case == "ic-block-drop") {
    double length = ic_vm["length"].as<double>();
    double width = ic_vm["width"].as<double>();
    // Error handling for the block size (length, width)
    if (length <= 0 || width <= 0) {
      std::cerr << "Error: Length and width must be positive!" << std::endl;
      exit(1);
    }
    double center_x = ic_vm["center_x"].as<double>();
    double center_y = ic_vm["center_y"].as<double>();
    // Error handling for the block initial position (center_x, center_y)
    if (center_x - length / 2.0 < left_wall ||
        center_x + length / 2.0 > right_wall ||
        center_y - width / 2.0 < bottom_wall ||
        center_y + width / 2.0 > top_wall) {
      std::cerr
          << "Error: The block must be within the domain boundaries! Please "
          << "adjust the center coordinates." << std::endl;
      exit(1);
    }
    sph = ic_block_drop(nb_particles, length, width, center_x, center_y);
  } else if (ic_case == "ic-droplet") {
    double radius = ic_vm["radius"].as<double>();
    // Error handling for the droplet radius
    if (radius <= 0) {
      std::cerr << "Error: Radius must be positive!" << std::endl;
      exit(1);
    }
    double center_x = ic_vm["center_x"].as<double>();
    double center_y = ic_vm["center_y"].as<double>();
    // Error handling for the droplet initial position (center_x, center_y)
    if (center_x - radius < left_wall || center_x + radius > right_wall ||
        center_y - radius < bottom_wall || center_y + radius > top_wall) {
      std::cerr << "Error: The droplet must be within the domain boundaries! "
                << "Please adjust the center coordinates." << std::endl;
      exit(1);
    }
    sph = ic_droplet(nb_particles, radius, center_x, center_y);
  } else {
    std::cerr << "Error: Function not found!" << std::endl;
  }

  po::variables_map constants_vm;
  std::ifstream constantsFile;
  constantsFile.open("../inputs/constants.txt");

  if (constantsFile.is_open()) {
    po::store(po::parse_config_file(constantsFile, desc), constants_vm);
    constantsFile.close();
  } else {
    std::cerr << "Error opening file: constants.txt" << std::endl;
  }

  po::notify(constants_vm);

  h = constants_vm["h"].as<double>();  // Radius of influence
  gas_constant = constants_vm["gas_constant"].as<double>();
  density_resting = constants_vm["density_resting"].as<double>();
  viscosity = constants_vm["viscosity"].as<double>();
  acceleration_gravity = constants_vm["acceleration_gravity"].as<double>();
  coeff_restitution = constants_vm["coeff_restitution"].as<double>();

  sph.set_timestep(dt);

  sph.set_rad_infl(h);
  sph.set_gas_constant(gas_constant);
  sph.set_density_resting(density_resting);
  sph.set_viscosity(viscosity);
  sph.set_acceleration_gravity(acceleration_gravity);
  sph.set_coeff_restitution(coeff_restitution);

  sph.set_left_wall(left_wall);
  sph.set_right_wall(right_wall);
  sph.set_bottom_wall(bottom_wall);
  sph.set_top_wall(top_wall);

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
