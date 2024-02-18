// This is a main program to run the SPH simulation
#include <filesystem>
#include <iostream>

#include "initial_conditions.h"
#include "main_prog_funcs.h"

// Start of the main programme
int main(int argc, char* argv[]) {
  // Define the output folder
  std::string OUTPUT_FOLDER = "../output";

  // Initialise an SphSolver object
  SphSolver sphSolver;

  // Declaring a fluid pointer (This is a nullptr at this stage)
  std::unique_ptr<Fluid> sphFluid;

  // Call the initialise function to initialise the objects based on the inputs
  initialise(sphFluid, sphSolver);

  std ::cout << "Initialisation finished -- OK"
             << "\n";

  // Initialise output files
  auto [initialPositions, finalPositionsFile, energiesFile] =
      initOutputFiles(OUTPUT_FOLDER);

  std ::cout << "Output files created -- OK"
             << "\n";

  // Store particles' positions before integration has started
  storeToFile(*sphFluid, "position", initialPositions);

  // Time integration loop
  sphSolver.timeIntegration(*sphFluid, finalPositionsFile, energiesFile);

  std ::cout << "SPH-SOLVER executed successfully -- OK"
             << "\n";

  return 0;
}

void initialise(std::unique_ptr<Fluid>& fluidPtr, SphSolver& sphSolver) {
  // Process to obtain the inputs provided by the user
  po::options_description desc("Allowed options");
  desc.add_options()("init_condition", po::value<std::string>(),
                     "take an initial condition")("T", po::value<double>(),
                                                  "take integration time")(
      "dt", po::value<double>(), "take time-step")(
      "coeffCfl1", po::value<double>(), "take lamda nu")(
      "coeffCfl2", po::value<double>(), "take lamda f")(
      "adaptive_timestep", po::value<bool>(),
      "take flag for adaptive time-step")("h", po::value<double>(),
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

  // Map the inputs read from the input files to expected inputs
  std::string icCase;
  std::vector<std::string> fileNames = {"case.txt", "domain.txt",
                                        "constants.txt"};
  po::variables_map caseVm, domainVm, constantsVm, icVm;
  std::vector<po::variables_map*> variableMaps = {&caseVm, &domainVm,
                                                  &constantsVm};

  for (size_t i = 0; i < fileNames.size(); ++i) {
    retrieveInputsFromFile(fileNames[i], icCase, desc, *(variableMaps[i]));
  }

  icCase = caseVm["init_condition"].as<std::string>();
  std::string fileName = icCase + ".txt";
  retrieveInputsFromFile(fileName, icCase, desc, icVm);

  handleInputErrors(caseVm, domainVm, constantsVm, icVm);

  setInitialConditions(icCase, fluidPtr, icVm, domainVm);

  // Set the parameters of the solver for the specific simulation
  sphSolver.setTimestep(caseVm["dt"].as<double>());
  sphSolver.setTotalIterations(
      ceil(caseVm["T"].as<double>() / caseVm["dt"].as<double>()));
  sphSolver.setOutputFrequency(caseVm["output_frequency"].as<int>());
  sphSolver.setCoeffRestitution(constantsVm["coeff_restitution"].as<double>());
  sphSolver.setLeftWall(domainVm["left_wall"].as<double>());
  sphSolver.setRightWall(domainVm["right_wall"].as<double>());
  sphSolver.setTopWall(domainVm["top_wall"].as<double>());
  sphSolver.setBottomWall(domainVm["bottom_wall"].as<double>());
  sphSolver.setPrecalculatedValues(constantsVm["h"].as<double>());

  // Define the fluid based on the inputs
  fluidPtr->setRadInfl(constantsVm["h"].as<double>());
  fluidPtr->setGasConstant(constantsVm["gas_constant"].as<double>());
  fluidPtr->setDensityResting(constantsVm["density_resting"].as<double>());
  fluidPtr->setViscosity(constantsVm["viscosity"].as<double>());
  fluidPtr->setAccelerationGravity(
      constantsVm["acceleration_gravity"].as<double>());

  // Calculate the mass of the particles
  sphSolver.createGrid(*fluidPtr);
  sphSolver.neighbourParticlesSearch(*fluidPtr);
  std::vector<std::vector<std::pair<int, double>>> neighboursPerParticle =
      sphSolver.getNeighbourParticles();
  fluidPtr->calculateMass(neighboursPerParticle);

  return;
}

void retrieveInputsFromFile(std::string fileName, std::string icCase,
                            po::options_description desc,
                            po::variables_map& vm) {
  std::ifstream caseFile;
  std::string errorMessage = "Error opening file: " + fileName;
  if (fileName == icCase + ".txt") {
    errorMessage +=
        " Make sure that the value of the init_condition in the case.txt "
        "file is one of the following: ic-one-particle, ic-two-particles, "
        "ic-three-particles, ic-four-particles, ic-droplet, ic-block-drop.";
  }
  // Try to open the file
  try {
    caseFile.open("../input/" + fileName);
    // Throw an exception if the file cannot be opened
    if (!caseFile.is_open()) {
      throw std::runtime_error(errorMessage);
    }
    po::store(po::parse_config_file(caseFile, desc), vm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(vm);
}

void handleInputErrors(po::variables_map caseVm, po::variables_map domainVm,
                       po::variables_map constantsVm, po::variables_map icVm) {
  try {
    // Error handling for the total integration time
    if (caseVm["T"].as<double>() <= 0) {
      throw std::runtime_error(
          "Error: Total integration time must be positive!");
      // Error handling for the time step
    } else if (caseVm["dt"].as<double>() <= 0 or
               caseVm["dt"].as<double>() > caseVm["T"].as<double>()) {
      throw std::runtime_error(
          "Error: Time step must be positive and lower than the total "
          "integration time!");
      // Error handling for the output frequency
    } else if (caseVm["output_frequency"].as<int>() <= 0 or
               caseVm["output_frequency"].as<int>() >
                   ceil(caseVm["T"].as<double>() / caseVm["dt"].as<double>())) {
      throw std::runtime_error(
          "Error: Output frequency must be positive and lower than the total "
          "number of iterations!");
      // Error handling for the domain boundaries input
    else if (caseVm["coeffCfl1"].as<double>() <= 0 or
        caseVm["coeffCfl1"].as<double>() >= 1 or
        caseVm["coeffCfl2"].as<double>() <= 0 or
        caseVm["coeffCfl2"].as<double>() >= 1) {
      throw std::runtime_error(
          "Error: The CFL coefficients must be positive and less than 1");
    } else if (domainVm["left_wall"].as<double>() >=
                   domainVm["right_wall"].as<double>() ||
               domainVm["bottom_wall"].as<double>() >=
                   domainVm["top_wall"].as<double>()) {
      throw std::runtime_error(
          "Error: Please adjust your domain boundaries so that left_wall < "
          "right wall and bottom_wall < top_wall.");
      // Error handling for the number of particles
    } else if (icVm["n"].as<int>() <= 0) {
      throw std::runtime_error("Error: Number of particles must be positive!");
    }
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
}

void setInitialConditions(std::string icCase, std::unique_ptr<Fluid>& fluidPtr,
                          po::variables_map icVm, po::variables_map domainVm) {
  // Fixed nbParticles ic cases map
  std::map<std::string, int> initConditionToParticlesMap = {
      {"ic-one-particle", 1},
      {"ic-two-particles", 2},
      {"ic-three-particles", 3},
      {"ic-four-particles", 4}};

  // Get the number of particles based on the ic case (for the more complex ic)
  int nbParticles;
  if (icCase == "ic-droplet" || icCase == "ic-block-drop") {
    nbParticles = icVm["n"].as<int>();
  } else {
    nbParticles = initConditionToParticlesMap[icCase];
  }

  // Fixed particles ic cases
  if (icCase == "ic-one-particle" || icCase == "ic-two-particles" ||
      icCase == "ic-three-particles" || icCase == "ic-four-particles") {
    std::vector<double> initX;
    initX.reserve(nbParticles);
    std::vector<double> initY;
    initY.reserve(nbParticles);

    for (int i = 0; i < nbParticles; i++) {
      initX[i] = icVm["init_x_" + std::to_string(i + 1)].as<double>();
      initY[i] = icVm["init_y_" + std::to_string(i + 1)].as<double>();
      // Error handling for the initial positions
      try {
        if (initX[i] < domainVm["left_wall"].as<double>() ||
            initX[i] > domainVm["right_wall"].as<double>() ||
            initY[i] < domainVm["bottom_wall"].as<double>() ||
            initY[i] > domainVm["top_wall"].as<double>()) {
          throw std::runtime_error(
              "Error: Particles must be within the domain boundaries! Please "
              "adjust the initial position coordinates.");
        }
      } catch (std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
      }
    }
    icBasic(fluidPtr, nbParticles, initX, initY);

    // Block drop case
  } else if (icCase == "ic-block-drop") {
    // Get the block dimensions and center coordinates from the ic file
    double length = icVm["length"].as<double>();
    double width = icVm["width"].as<double>();
    // Error handling for the block size (length, width)
    try {
      if (length <= 0 || width <= 0) {
        throw std::runtime_error("Error: Length and width must be positive!");
      }
    } catch (std::runtime_error& e) {
      std::cerr << e.what() << std::endl;
      exit(1);
    }
    double centerX = icVm["center_x"].as<double>();
    double centerY = icVm["center_y"].as<double>();
    // Error handling for the block initial position (centerX, centerY)
    try {
      if (centerX - length / 2.0 < domainVm["left_wall"].as<double>() ||
          centerX + length / 2.0 > domainVm["right_wall"].as<double>() ||
          centerY - width / 2.0 < domainVm["bottom_wall"].as<double>() ||
          centerY + width / 2.0 > domainVm["top_wall"].as<double>()) {
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
    icBlockDrop(fluidPtr, nbParticles, length, width, centerX, centerY);
    // Droplet case
  } else if (icCase == "ic-droplet") {
    // Get the droplet radius and center coordinates from the ic file
    double radius = icVm["radius"].as<double>();
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
    double centerX = icVm["center_x"].as<double>();
    double centerY = icVm["center_y"].as<double>();
    // Error handling for the droplet initial position (centerX, centerY)
    try {
      if (centerX - radius < domainVm["left_wall"].as<double>() ||
          centerX + radius > domainVm["right_wall"].as<double>() ||
          centerY - radius < domainVm["bottom_wall"].as<double>() ||
          centerY + radius > domainVm["top_wall"].as<double>()) {
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
    icDroplet(fluidPtr, nbParticles, radius, centerX, centerY);
  } else {
    std::cerr << "Error: Initial condition function not found! Make sure "
              << "that the value of the init_condition in the case.txt file is "
              << "one of the following: ic-one-particle, ic-two-particles, "
              << "ic-three-particles, ic-four-particles, ic-droplet, "
              << "ic-block-drop." << std::endl;
    exit(1);
  }

  // Map the inputs read from the constants file to expected inputs
  po::variables_map constantsVm;
  std::ifstream constantsFile;
  try {
    constantsFile.open("../input/constants.txt");
    // Throw an exception if the file cannot be opened
    if (!constantsFile.is_open()) {
      throw std::runtime_error("Error opening file: constants.txt");
    }
    po::store(po::parse_config_file(constantsFile, desc), constantsVm);
  } catch (std::runtime_error& e) {
    // Handle the exception by printing the error message and exiting the
    // program
    std::cerr << e.what() << std::endl;
    exit(1);
  }
  po::notify(constantsVm);

  // Set the parameters of the solver for the specific simulation
  sphSolver.setTimestep(caseVm["dt"].as<double>());
  sphSolver.setAdaptiveTimestep(caseVm["adaptive_timestep"].as<bool>());
  sphSolver.setCflCoefficients(caseVm["coeffCfl1"].as<double>(),
                               caseVm["coeffCfl2"].as<double>());
  sphSolver.setTotalTime(totalTime);
  sphSolver.setOutputFrequency(caseVm["output_frequency"].as<int>());
  sphSolver.setCoeffRestitution(constantsVm["coeff_restitution"].as<double>());
  sphSolver.setLeftWall(domainVm["left_wall"].as<double>());
  sphSolver.setRightWall(domainVm["right_wall"].as<double>());
  sphSolver.setTopWall(domainVm["top_wall"].as<double>());
  sphSolver.setBottomWall(domainVm["bottom_wall"].as<double>());
  sphSolver.setPrecalculatedValues(constantsVm["h"].as<double>());

  // Define the fluid based on the inputs
  fluidPtr->setRadInfl(constantsVm["h"].as<double>());
  fluidPtr->setGasConstant(constantsVm["gas_constant"].as<double>());
  fluidPtr->setDensityResting(constantsVm["density_resting"].as<double>());
  fluidPtr->setViscosity(constantsVm["viscosity"].as<double>());
  fluidPtr->setAccelerationGravity(
      constantsVm["acceleration_gravity"].as<double>());

  // Calculate the mass of the particles
  sphSolver.createGrid(*fluidPtr);
  sphSolver.neighbourParticlesSearch(*fluidPtr);
  std::vector<std::vector<std::pair<int, double>>> neighboursPerParticle =
      sphSolver.getNeighbourParticles();
  fluidPtr->calculateMass(neighboursPerParticle);

  return;
}

void createDirectory(std::string folderPath) {
  // Check if the target folder already exists
  if (!std::filesystem::exists(folderPath)) {
    // Create target folder
    std::filesystem::create_directories(folderPath);
  }
}

std::tuple<std::ofstream, std::ofstream, std::ofstream> initOutputFiles(
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
  energies << "dt,t,Ek,Ep,Etotal"
           << "\n";

  return std::make_tuple(std::move(initialPositions), std::move(finalPositions),
                         std::move(energies));
}

void storeToFile(Fluid& fluid, std::string type, std::ofstream& targetFile,
                 double dt, double currentTime) {
  if (type == "energy") {
    // Write energies in the Energy-File
    targetFile << dt << "," << currentTime << "," << fluid.getKineticEnergy()
               << "," << fluid.getPotentialEnergy() << ","
               << fluid.getPotentialEnergy() + fluid.getKineticEnergy() << "\n";
  } else if (type == "position") {
    // Write positions in the position file
    for (int k = 0; k < fluid.getNumberOfParticles(); k++) {
      targetFile << fluid.getPositionX(k) << "," << fluid.getPositionY(k)
                 << "\n";
    }
  }
}
