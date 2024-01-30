// This is a main program to run the SPH simulation
#include <boost/program_options.hpp>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>

#include "initial_conditions.h"
#include "main_prog_funcs.h"

// Start of the main programme
int main(int argc, char* argv[]) {
  std::string OUTPUT_FOLDER = "../output";
  // Read input files, initialise the sph class and the parameters of the
  // problem
  sph_solver sph_solver;

  fluid* sph_fluid = new fluid;  // Allocate memory for the fluid object

  initialise(&sph_fluid, sph_solver);

  std ::cout << "Initialisation finished -- OK"
             << "\n";

  // Initialise output files
  auto [initialPositions, finalPositionsFile, energiesFile] =
      init_output_files(OUTPUT_FOLDER);

  std ::cout << "Output files created -- OK"
             << "\n";

  // Store particles' positions before integration has started
  storeToFile(*sph_fluid, "position", initialPositions);

  // Time integration loop
  sph_solver.time_integration(*sph_fluid, finalPositionsFile, energiesFile);

  std ::cout << "SPH-SOLVER executed successfully -- OK"
             << "\n";

  return 0;
}

void initialise(fluid** fluid_ptr, sph_solver& sph_solver) {
  int nb_particles;  //  Number of particles

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
      "output_frequency", po::value<int>(),
      "take frequency that output will be written to file");

  // Map the inputs read from the case file to expected inputs
  po::variables_map case_vm;
  std::ifstream caseFile;

  // Try to open the case.txt file
  try {
    caseFile.open("../input/case.txt");
    // Throw an exception if the file cannot be opened
    if (!caseFile.is_open()) {
      throw std::runtime_error("Error opening file: case.txt");
    }
    po::store(po::parse_config_file(caseFile, desc), case_vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(case_vm);

  double total_time = case_vm["T"].as<double>();  // Total integration time
  // Error handling for the total integration time
  try {
    if (total_time <= 0) {
      throw std::runtime_error(
          "Error: Total integration time must be positive!");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  // Time step dt
  // Error handling for the time step
  try {
    if (case_vm["dt"].as<double>() <= 0 or
        case_vm["dt"].as<double>() > total_time) {
      throw std::runtime_error(
          "Error: Time step must be positive and lower than the total "
          "integration time!");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  // Error handling for the output frequency
  try {
    if (case_vm["output_frequency"].as<int>() <= 0 or
        case_vm["output_frequency"].as<int>() >
            ceil(total_time / case_vm["dt"].as<double>())) {
      throw std::runtime_error(
          "Error: Output frequency must be positive and lower than the total "
          "number of iterations!");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  // Map the inputs read from the domain file to expected inputs
  po::variables_map domain_vm;
  std::ifstream domainFile;
  try {
    domainFile.open("../input/domain.txt");
    // Throw an exception if the file cannot be opened
    if (!domainFile.is_open()) {
      throw std::runtime_error("Error opening file: domain.txt");
    }
    po::store(po::parse_config_file(domainFile, desc), domain_vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(domain_vm);

  // Error handling for the domain boundaries input
  try {
    if (domain_vm["left_wall"].as<double>() >=
            domain_vm["right_wall"].as<double>() ||
        domain_vm["bottom_wall"].as<double>() >=
            domain_vm["top_wall"].as<double>()) {
      throw std::runtime_error(
          "Error: Please adjust your domain boundaries so that left_wall < "
          "right wall and bottom_wall < top_wall.");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }

  // Map the inputs read from the initial condition file to expected inputs
  std::string ic_case = case_vm["init_condition"].as<std::string>();
  po::variables_map ic_vm;
  std::ifstream icFile;
  // Open the file of the initial condition the user has chosen
  try {
    icFile.open("../input/" + ic_case + ".txt");
    // Throw an exception if the file cannot be opened
    if (!icFile.is_open()) {
      throw std::runtime_error(
          "Error opening file: " + ic_case +
          ".txt Make sure that the value of the init_condition in the case.txt "
          "file is one of the following: ic-one-particle, ic-two-particles, "
          "ic-three-particles, ic-four-particles, ic-droplet, ic-block-drop.");
    }
    po::store(po::parse_config_file(icFile, desc), ic_vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
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
    try {
      if (nb_particles <= 0) {
        throw std::runtime_error(
            "Error: Number of particles must be positive!");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
  } else {
    nb_particles = initConditionToParticlesMap[case_vm["init_condition"]
                                                   .as<std::string>()];
  }

  // Fixed particles ic cases
  if (ic_case == "ic-one-particle" || ic_case == "ic-two-particles" ||
      ic_case == "ic-three-particles" || ic_case == "ic-four-particles") {
    // Get particles' initial poistions from the ic file
    double* init_x = new double[nb_particles];
    double* init_y = new double[nb_particles];
    for (int i = 0; i < nb_particles; i++) {
      init_x[i] = ic_vm["init_x_" + std::to_string(i + 1)].as<double>();
      init_y[i] = ic_vm["init_y_" + std::to_string(i + 1)].as<double>();
      // Error handling for the initial positions
      try {
        if (init_x[i] < domain_vm["left_wall"].as<double>() ||
            init_x[i] > domain_vm["right_wall"].as<double>() ||
            init_y[i] < domain_vm["bottom_wall"].as<double>() ||
            init_y[i] > domain_vm["top_wall"].as<double>()) {
          throw std::runtime_error(
              "Error: Particles must be within the domain boundaries! Please "
              "adjust the initial position coordinates.");
        }
      } catch (std::runtime_error& e) {
        // Handle the exception by printing the error message and exiting the
        // program
        std::cerr << e.what() << std::endl;
        delete[] init_x;
        delete[] init_y;
        exit(1);
      }
    }
    ic_basic(fluid_ptr, nb_particles, init_x, init_y);
    delete[] init_x;
    delete[] init_y;
    // Block drop case
  } else if (ic_case == "ic-block-drop") {
    // Get the block dimensions and center coordinates from the ic file
    double length = ic_vm["length"].as<double>();
    double width = ic_vm["width"].as<double>();
    // Error handling for the block size (length, width)
    try {
      if (length <= 0 || width <= 0) {
        throw std::runtime_error("Error: Length and width must be positive!");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    double center_x = ic_vm["center_x"].as<double>();
    double center_y = ic_vm["center_y"].as<double>();
    // Error handling for the block initial position (center_x, center_y)
    try {
      if (center_x - length / 2.0 < domain_vm["left_wall"].as<double>() ||
          center_x + length / 2.0 > domain_vm["right_wall"].as<double>() ||
          center_y - width / 2.0 < domain_vm["bottom_wall"].as<double>() ||
          center_y + width / 2.0 > domain_vm["top_wall"].as<double>()) {
        throw std::runtime_error(
            "Error: The block must be within the domain boundaries! Please "
            "adjust the center coordinates.");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    ic_block_drop(fluid_ptr, nb_particles, length, width, center_x, center_y);
    // Droplet case
  } else if (ic_case == "ic-droplet") {
    // Get the droplet radius and center coordinates from the ic file
    double radius = ic_vm["radius"].as<double>();
    // Error handling for the droplet radius
    try {
      if (radius <= 0) {
        throw std::runtime_error("Error: Radius must be positive!");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    double center_x = ic_vm["center_x"].as<double>();
    double center_y = ic_vm["center_y"].as<double>();
    // Error handling for the droplet initial position (center_x, center_y)
    try {
      if (center_x - radius < domain_vm["left_wall"].as<double>() ||
          center_x + radius > domain_vm["right_wall"].as<double>() ||
          center_y - radius < domain_vm["bottom_wall"].as<double>() ||
          center_y + radius > domain_vm["top_wall"].as<double>()) {
        throw std::runtime_error(
            "Error: The droplet must be within the domain boundaries! Please "
            "adjust the center coordinates.");
      }
    } catch (std::runtime_error& e) {
      // Handle the exception by printing the error message and exiting the
      // program
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    ic_droplet(fluid_ptr, nb_particles, radius, center_x, center_y);
  } else {
    std::cerr << "Error: Initial condition function not found! Make sure "
              << "that the value of the init_condition in the case.txt file is "
              << "one of the following: ic-one-particle, ic-two-particles, "
              << "ic-three-particles, ic-four-particles, ic-droplet, "
              << "ic-block-drop." << std::endl;
    exit(1);
  }

  // Map the inputs read from the constants file to expected inputs
  po::variables_map constants_vm;
  std::ifstream constantsFile;
  try {
    constantsFile.open("../input/constants.txt");
    // Throw an exception if the file cannot be opened
    if (!constantsFile.is_open()) {
      throw std::runtime_error("Error opening file: constants.txt");
    }
    po::store(po::parse_config_file(constantsFile, desc), constants_vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(constants_vm);

  sph_solver.set_timestep(case_vm["dt"].as<double>());
  sph_solver.set_total_iter(ceil(total_time / case_vm["dt"].as<double>()));
  sph_solver.set_output_frequency(case_vm["output_frequency"].as<int>());
  sph_solver.set_coeff_restitution(
      constants_vm["coeff_restitution"].as<double>());
  sph_solver.set_left_wall(domain_vm["left_wall"].as<double>());
  sph_solver.set_right_wall(domain_vm["right_wall"].as<double>());
  sph_solver.set_top_wall(domain_vm["top_wall"].as<double>());
  sph_solver.set_bottom_wall(domain_vm["bottom_wall"].as<double>());

  fluid* objPtr = *fluid_ptr;

  objPtr->set_rad_infl(constants_vm["h"].as<double>());
  objPtr->set_gas_constant(constants_vm["gas_constant"].as<double>());
  objPtr->set_density_resting(constants_vm["density_resting"].as<double>());
  objPtr->set_viscosity(constants_vm["viscosity"].as<double>());
  objPtr->set_acceleration_gravity(
      constants_vm["acceleration_gravity"].as<double>());

  // Calculate the mass of the particles
  objPtr->calc_mass();

  return;
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

void storeToFile(fluid& fluid, std::string type, std::ofstream& targetFile,
                 double dt, int currentIteration) {
  if (type == "energy") {
    // Write energies on the Energy-File
    targetFile << currentIteration * dt << "," << fluid.get_kinetic_energy()
               << "," << fluid.get_potential_energy() << ","
               << fluid.get_potential_energy() + fluid.get_kinetic_energy()
               << "\n";
  } else if (type == "position") {
    for (int k = 0; k < fluid.get_number_of_particles(); k++) {
      targetFile << fluid.get_position_x(k) << "," << fluid.get_position_y(k)
                 << "\n";
    }
  }
}
